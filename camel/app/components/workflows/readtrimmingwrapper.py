from dataclasses import dataclass
from typing import List, Dict, Any

import os

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import SNAKEFILE_READ_TRIMMING
from camel.resources.snakefile.read_trimming import OUTPUT_READ_TRIMMING_REPORT, OUTPUT_READ_TRIMMING_READS_PE, \
    OUTPUT_READ_TRIMMING_READS_SE_REV, OUTPUT_READ_TRIMMING_READS_SE_FWD, OUTPUT_READ_TRIMMING_INFORMS, \
    OUTPUT_READ_TRIMMING_FASTQC_PRE


class ReadTrimmingWrapper(object):
    """
    This class is used as a wrapper class around the read trimming Snakemake workflow.
    """

    @dataclass
    class ReadTrimmingOutput:
        report_section: HtmlReportSection
        trimmed_reads_pe: List[ToolIOFile]
        trimmed_reads_se_fwd: List[ToolIOFile]
        trimmed_reads_se_rev: List[ToolIOFile]
        informs_trimmomatic: Dict[str, Any]
        fastq_reports_pre: List[ToolIOFile]

    def __init__(self, working_dir: str) -> None:
        """
        Initializes the read trimming wrapper.
        :param working_dir: Working directory
        """
        self._working_dir = working_dir
        self._output = None

    def run_workflow(self, pe_reads: List[str], threads: int = 8, export_fastq: bool = False) -> None:
        """
        Runs the read trimming workflow.
        :param pe_reads: Input PE FASTQ reads
        :param threads: Number of threads to use
        :param export_fastq: If True, FASTQ files are included in the report
        :return: None
        """
        config_data = {
            'working_dir': self._working_dir,
            'fastq_pe': [
                {'name': os.path.basename(pe_reads[0]), 'path': pe_reads[0]},
                {'name': os.path.basename(pe_reads[1]), 'path': pe_reads[1]}],
            'read_trimming': {'export_fastq': str(export_fastq)}
        }
        output_path = os.path.join(self._working_dir, OUTPUT_READ_TRIMMING_REPORT)
        SnakePipelineUtils.run_snakemake(
            SNAKEFILE_READ_TRIMMING, config_data, [output_path], self._working_dir, threads)
        self.__set_output(output_path)

    def __set_output(self, output_path: str) -> None:
        """
        Sets the output of this tool.
        :param output_path: Report output path
        :return: None
        """
        self._output = ReadTrimmingWrapper.ReadTrimmingOutput(
            report_section=SnakemakeUtils.load_object(output_path)[0].value,
            trimmed_reads_pe=SnakemakeUtils.load_object(
                os.path.join(self._working_dir, OUTPUT_READ_TRIMMING_READS_PE)),
            trimmed_reads_se_fwd=SnakemakeUtils.load_object(
                os.path.join(self._working_dir, OUTPUT_READ_TRIMMING_READS_SE_FWD)),
            trimmed_reads_se_rev=SnakemakeUtils.load_object(
                os.path.join(self._working_dir, OUTPUT_READ_TRIMMING_READS_SE_REV)),
            informs_trimmomatic=SnakemakeUtils.load_object(
                os.path.join(self._working_dir, OUTPUT_READ_TRIMMING_INFORMS)),
            fastq_reports_pre=SnakemakeUtils.load_object(
                os.path.join(self._working_dir, OUTPUT_READ_TRIMMING_FASTQC_PRE))
        )

    @property
    def output(self) -> 'ReadTrimmingWrapper.ReadTrimmingOutput':
        """
        Returns the report section generated during the read trimming pipeline.
        :return: Trimming section
        """
        return self._output
