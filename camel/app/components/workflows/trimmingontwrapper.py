from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Optional

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import trimming_ont
from camel.resources.snakefile.trimming_ont import INPUT_ONT_FASTQ


class TrimmingONTWrapper(object):
    """
    This class is used as a wrapper class around the ONT read trimming Snakemake workflow.
    """

    @dataclass
    class ReadTrimmingOutput:
        report_section: HtmlReportSection
        tsv_summary: Path
        trimmed_reads: List[ToolIOFile]
        log_file: Optional[Path] = None

    def __init__(self, working_dir: Path) -> None:
        """
        Initializes the read trimming wrapper.
        :param working_dir: Working directory
        """
        self._working_dir = Path(working_dir)
        self._output = None

    def run_workflow(self, se_reads: Path, export_fastq: bool = False, threads: int = 8) -> None:
        """
        Runs the read trimming workflow.
        :param se_reads: Input SE FASTQ reads
        :param threads: Number of threads to use
        :param export_fastq: If True, FASTQ files are included in the report
        :return: None
        """
        path_io = self._working_dir / INPUT_ONT_FASTQ
        path_io.parent.mkdir(exist_ok=True, parents=True)
        SnakemakeUtils.dump_object([ToolIOFile(se_reads)], path_io)
        config_data = {
            'working_dir': str(self._working_dir),
            'read_trimming': {'export_fastq': str(export_fastq)},
            'read_type': 'ont'
        }
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._working_dir)
        output_files = {
            'HTML': self._working_dir / trimming_ont.OUTPUT_TRIMMING_ONT_REPORT,
            'TSV': self._working_dir / trimming_ont.OUTPUT_TRIMMING_ONT_SUMMARY,
        }
        SnakePipelineUtils.run_snakemake(
            trimming_ont.SNAKEFILE_TRIMMING_ONT, config_file, list(output_files.values()),
            self._working_dir, threads)
        self.__set_output(output_files)

    def __set_output(self, output_files: Dict[str, Path]) -> None:
        """
        Sets the output of this tool.
        :param output_files: Output files by key.
        :return: None
        """
        log_path = self._working_dir / 'camel.log'
        self._output = TrimmingONTWrapper.ReadTrimmingOutput(
            report_section=SnakemakeUtils.load_object(output_files['HTML'])[0].value,
            tsv_summary=output_files['TSV'],
            trimmed_reads=SnakemakeUtils.load_object(
                self._working_dir / trimming_ont.OUTPUT_TRIMMING_ONT_READS),
            log_file=log_path if log_path.exists() else None
        )

    @property
    def output(self) -> 'TrimmingONTWrapper.ReadTrimmingOutput':
        """
        Returns the report section generated during the read trimming pipeline.
        :return: Trimming section
        """
        return self._output
