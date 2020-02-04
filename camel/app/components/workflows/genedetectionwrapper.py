from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Any, List, Optional

from camel.app.components.genedetection.genedetectionhitbase import GeneDetectionHitBase
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import gene_detection


@dataclass
class GeneDetectionOutput:
    report_section: HtmlReportSection
    detected_hits: List[GeneDetectionHitBase]
    informs: Dict[str, Any]
    log_file: Optional[Path] = None


class GeneDetectionWrapper(object):
    """
    This class is used as a wrapper class around the gene detection Snakemake workflow.
    """

    def __init__(self, working_dir: str) -> None:
        """
        Initializes the read trimming helper.
        :param working_dir: Working directory
        """
        self._working_dir = Path(working_dir)
        self._output = None

    def __run_workflow(self, config_data: Dict[str, Any], threads: int) -> None:
        """
        Runs the gene detection workflow with the given config data and number of threads.
        :param config_data: Config data
        :param threads: Number of threads
        :return: None
        """
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._working_dir)
        output_files = {
            'report': self._working_dir / str(gene_detection.OUTPUT_GENE_DETECTION_REPORT).format(db='db'),
            'informs': self._working_dir / str(gene_detection.OUTPUT_GENE_DETECTION_INFORMS).format(db='db'),
            'hits': self._working_dir / str(gene_detection.OUTPUT_GENE_DETECTION_ALL_HITS).format(db='db')
        }
        SnakePipelineUtils.run_snakemake(
            gene_detection.SNAKEFILE_GENE_DETECTION, config_file, list(output_files.values()), self._working_dir,
            threads)
        self.__set_output(output_files)

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
        self.__run_workflow(config_data, threads)

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
        self.__run_workflow(config_data, threads)

    def run_workflow_kma(self, fastq_pe_path: List[str], sample_name: str, db_data: Dict[str, Any],
                         threads: int = 8) -> None:
        """
        Runs the gene detection workflow using KMA.
        :param fastq_pe_path: Input PE FASTQ files
        :param sample_name: Sample name
        :param db_data: Database configuration
        :param threads: Number of threads to use
        :return: None
        """
        self.__create_input_srst2(fastq_pe_path)
        config_data = self.__get_config_data(sample_name, db_data, 'kma')
        self.__run_workflow(config_data, threads)

    def __create_input_blast(self, fasta_path: str) -> None:
        """
        Creates the input for the workflow in 'blast' mode.
        :param fasta_path: FASTA file path
        :return: None
        """
        path = self._working_dir / str(gene_detection.INPUT_GENE_DETECTION_FASTA).format(db='db')
        if not path.parent.exists():
            path.parent.mkdir(parents=True)
        SnakemakeUtils.dump_object([ToolIOFile(fasta_path)], str(path))

    def __create_input_srst2(self, fastq_pe: List[str]) -> None:
        """
        Creates the input for the workflow in 'srst2' mode.
        :param fastq_pe: Pair of FASTQ PE files
        :return: None
        """
        path = self._working_dir / 'fq_dict.io'
        fq_dict = {'PE': [ToolIOFile(x) for x in fastq_pe]}
        SnakemakeUtils.dump_object(fq_dict, path)

    def __get_config_data(self, sample_name: str, db_data: Dict[str, Any], detection_method: str) -> Dict[str, Any]:
        """
        Returns the configuration data for Snakemake.
        :param sample_name: Sample name
        :param db_data: Database information
        :param detection_method: Detection method
        :return: Config data
        """
        return {
            'working_dir': str(self._working_dir),
            'sample_name': sample_name,
            'detection_method': detection_method,
            'gene_detection': {'db': db_data}
        }

    def __set_output(self, output_files: Dict[str, str]) -> None:
        """
        Sets the output of the workflow.
        :param output_files: Output files dictionary
        :return: None
        """
        log_file_path = self._working_dir / 'camel.log'
        self._output = GeneDetectionOutput(
            report_section=SnakemakeUtils.load_object(output_files['report'])[0].value,
            detected_hits=[v.value for v in SnakemakeUtils.load_object(output_files['hits'])],
            informs=SnakemakeUtils.load_object(output_files['informs']),
            log_file=log_file_path if log_file_path.exists() else None
        )

    @property
    def output(self) -> GeneDetectionOutput:
        """
        Returns the output generated by the gene detection workflow.
        :return: Output
        """
        return self._output
