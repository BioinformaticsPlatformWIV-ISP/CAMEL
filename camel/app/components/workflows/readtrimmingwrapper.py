from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Any, Optional

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import trimming_illumina


class ReadTrimmingWrapper(object):
    """
    This class is used as a wrapper class around the Illumina read trimming Snakemake workflow.
    """

    @dataclass
    class ReadTrimmingOutput:
        report_section: HtmlReportSection
        trimmed_reads_pe: List[ToolIOFile]
        trimmed_reads_se_fwd: List[ToolIOFile]
        trimmed_reads_se_rev: List[ToolIOFile]
        informs_trimmomatic: Dict[str, Any]
        fastq_reports_pre: List[ToolIOFile]
        log_file: Optional[Path] = None

    def __init__(self, working_dir: str) -> None:
        """
        Initializes the read trimming wrapper.
        :param working_dir: Working directory
        """
        self._working_dir = Path(working_dir)
        self._output = None

    def run_workflow(self, pe_reads: List[str], threads: int = 8, export_fastq: bool = False) -> None:
        """
        Runs the read trimming workflow.
        :param pe_reads: Input PE FASTQ reads
        :param threads: Number of threads to use
        :param export_fastq: If True, FASTQ files are included in the report
        :return: None
        """
        pe_reads_path = [Path(r) for r in pe_reads]
        config_data = {
            'working_dir': str(self._working_dir),
            'fastq_pe': [{'name': p.name, 'path': str(p)} for p in pe_reads_path],
            'read_trimming': {'export_fastq': str(export_fastq)}
        }
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._working_dir)
        output_path = self._working_dir / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_REPORT
        SnakePipelineUtils.run_snakemake(
            trimming_illumina.SNAKEFILE_TRIMMING_ILLUMINA, config_file, [output_path], self._working_dir, threads)
        self.__set_output(output_path)

    def __set_output(self, output_path: str) -> None:
        """
        Sets the output of this tool.
        :param output_path: Report output path
        :return: None
        """
        log_path = self._working_dir / 'camel.log'
        self._output = ReadTrimmingWrapper.ReadTrimmingOutput(
            report_section=SnakemakeUtils.load_object(output_path)[0].value,
            trimmed_reads_pe=SnakemakeUtils.load_object(str(
                self._working_dir / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_READS_PE)),
            trimmed_reads_se_fwd=SnakemakeUtils.load_object(str(
                self._working_dir / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_READS_SE_FWD)),
            trimmed_reads_se_rev=SnakemakeUtils.load_object(str(
                self._working_dir / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_READS_SE_REV)),
            informs_trimmomatic=SnakemakeUtils.load_object(str(
                self._working_dir / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_INFORMS)),
            fastq_reports_pre=SnakemakeUtils.load_object(str(
                self._working_dir / trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_FASTQC_PRE)),
            log_file=log_path if log_path.exists() else None
        )

    @property
    def output(self) -> 'ReadTrimmingWrapper.ReadTrimmingOutput':
        """
        Returns the report section generated during the read trimming pipeline.
        :return: Trimming section
        """
        return self._output
