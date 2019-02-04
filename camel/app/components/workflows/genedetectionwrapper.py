from dataclasses import dataclass
from typing import Dict, Any, List, Optional

import os

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import SNAKEFILE_GENE_DETECTION
from camel.resources.snakefile.gene_detection import INPUT_GENE_DETECTION_FASTA, OUTPUT_GENE_DETECTION_REPORT, \
    INPUT_GENE_DETECTION_FASTQ


class GeneDetectionWrapper(object):
    """
    This class is used as a wrapper class around the gene detection Snakemake workflow.
    """

    @dataclass
    class GeneDetectionOutput:
        report_section: HtmlReportSection
        log_file: Optional[str] = None

    def __init__(self, working_dir: str) -> None:
        """
        Initializes the read trimming helper.
        :param working_dir: Working directory
        """
        self._working_dir = working_dir
        self._output = None

    def run_workflow_blast(self, fasta_path: str, sample_name: str, db_data: Dict[str, Any], threads: int = 8) -> None:
        """
        Runs the gene detection workflow using BLAST.
        :param fasta_path: Input FASTA file
        :param sample_name: Sample name
        :param db_data: Database configuration data (should contain at least 'path' as key referring to the location of
          the database.)
        :param threads: Number of threads to use
        :return: None
        """
        self.__create_input_blast(fasta_path)
        config_data = self.__get_config_data(sample_name, db_data, 'blast')
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._working_dir)
        output_path = os.path.join(self._working_dir, OUTPUT_GENE_DETECTION_REPORT.format(db='db'))
        SnakePipelineUtils.run_snakemake(
            SNAKEFILE_GENE_DETECTION, config_file, [output_path], self._working_dir, threads)
        self.__set_output(output_path)

    def run_workflow_srst2(self, fastq_pe_path: List[str], sample_name: str, db_data: Dict[str, Any],
                           threads: int = 8) -> None:
        """
        Runs the gene detection workflow using SRST2.
        :param fastq_pe_path: Input PE FASTQ files
        :param sample_name: Sample name
        :param db_data: Database configuration
        :param threads: Number of threads to use
        :return: None
        """
        self.__create_input_srst2(fastq_pe_path)
        config_data = self.__get_config_data(sample_name, db_data, 'srst2')
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._working_dir)
        output_path = os.path.join(self._working_dir, OUTPUT_GENE_DETECTION_REPORT.format(db='db'))
        SnakePipelineUtils.run_snakemake(
            SNAKEFILE_GENE_DETECTION, config_file, [output_path], self._working_dir, threads)
        self.__set_output(output_path)

    def __create_input_blast(self, fasta_path: str) -> None:
        """
        Creates the input for the workflow in 'blast' mode.
        :param fasta_path: FASTA file path
        :return: None
        """
        path = os.path.join(self._working_dir, INPUT_GENE_DETECTION_FASTA.format(db='db'))
        target_dir = os.path.dirname(path)
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)
        SnakemakeUtils.dump_object([ToolIOFile(fasta_path)], path)

    def __create_input_srst2(self, fastq_pe: List[str]) -> None:
        """
        Creates the input for the workflow in 'srst2' mode.
        :param fastq_pe: Pair of FASTQ PE files
        :return: None
        """
        path = os.path.join(self._working_dir, INPUT_GENE_DETECTION_FASTQ.format(db='db'))
        target_dir = os.path.dirname(path)
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)
        SnakemakeUtils.dump_object([ToolIOFile(fastq_pe[0]), ToolIOFile(fastq_pe[1])], path)

    def __get_config_data(self, sample_name: str, db_data: Dict[str, Any], detection_method: str) -> Dict[str, Any]:
        """
        Returns the configuration data for Snakemake.
        :param sample_name: Sample name
        :param db_data: Database information
        :param detection_method: Detection method
        :return: Config data
        """
        return {
            'working_dir': self._working_dir,
            'sample_name': sample_name,
            'detection_method': detection_method,
            'gene_detection': {'db': db_data}
        }

    def __set_output(self, report_path: str) -> None:
        """
        Sets the output of the workflow.
        :param report_path: Report path
        :return: None
        """
        log_file_path = os.path.join(self._working_dir, 'camel.log')
        self._output = GeneDetectionWrapper.GeneDetectionOutput(
            report_section=SnakemakeUtils.load_object(report_path)[0].value,
            log_file=log_file_path if os.path.isfile(log_file_path) else None
        )

    @property
    def output(self) -> 'GeneDetectionWrapper.GeneDetectionOutput':
        """
        Returns the output generated by the gene detection workflow.
        :return: Output
        """
        return self._output
