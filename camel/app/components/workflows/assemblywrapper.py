from dataclasses import dataclass
from typing import List, Optional, Dict, Any, Union

import os

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import SNAKEFILE_ASSEMBLY_SPADES
from camel.resources.snakefile.assembly_spades import OUTPUT_ASSEMBLY_REPORT, OUTPUT_ASSEMBLY_FASTA, \
    OUTPUT_ASSEMBLY_INFORMS, OUTPUT_ASSEMBLY_FILTERING_INFORMS
from camel.resources.snakefile.read_trimming import OUTPUT_READ_TRIMMING_READS_PE, \
    OUTPUT_READ_TRIMMING_READS_SE_REV, OUTPUT_READ_TRIMMING_READS_SE_FWD


@dataclass
class AssemblyOutput:
    report_section: HtmlReportSection
    fasta_contigs: ToolIOFile
    log_file: Optional[str]
    informs: List[Dict[str, Any]]


class AssemblyWrapper(object):
    """
    This class is used as a wrapper class around the assembly Snakemake workflow.
    """

    def __init__(self, working_dir: str) -> None:
        """
        Initializes the read trimming helper.
        :param working_dir: Working directory
        """
        self._working_dir = working_dir
        self._output = None

    def run_workflow(self, sample_name: str, reads_pe: List[ToolIOFile], reads_se_fwd: List[ToolIOFile],
                     reads_se_rev: List[ToolIOFile], kmers: str = None, cov_cutoff: Union[str, int] = 'off',
                     min_contig_length: Optional[int] = None, threads: int = 8) -> None:
        """
        Runs the read trimming workflow.
        :param reads_pe: Paired end reads
        :param reads_se_fwd: Single end forward reads
        :param reads_se_rev: Single end reverse reads
        :param sample_name: Sample name
        :param kmers: Comma separated list of Kmer sizes to use for assembly
        :param cov_cutoff: Coverage cutoff
        :param min_contig_length: Minimum contig length
        :param threads: Number of threads
        :return: None
        """
        self.__prepare_input_files(reads_pe, reads_se_fwd, reads_se_rev)
        config_data = self.__get_config_data(sample_name, kmers, cov_cutoff, min_contig_length)
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._working_dir)
        output_files = {
            'HTML': os.path.join(self._working_dir, OUTPUT_ASSEMBLY_REPORT),
            'FASTA': os.path.join(self._working_dir, OUTPUT_ASSEMBLY_FASTA),
            'INFORMS_spades': os.path.join(self._working_dir, OUTPUT_ASSEMBLY_INFORMS),
        }
        if min_contig_length is not None:
            output_files['INFORMS_seqtk'] = os.path.join(self._working_dir, OUTPUT_ASSEMBLY_FILTERING_INFORMS)
        SnakePipelineUtils.run_snakemake(
            SNAKEFILE_ASSEMBLY_SPADES, config_file, list(output_files.values()), self._working_dir, threads)
        self.__set_output(output_files)

    def __get_config_data(self, sample_name: str, kmers: str, cov_cutoff: Union[int, str],
                          min_contig_length: Optional[int] = None) -> Dict[str, Any]:
        """
        Builds the configuration file to run the assembly workflow.
        :param sample_name: Sample name
        :param kmers: Comma separated list of Kmer sizes to use for the assembly
        :param cov_cutoff: Coverage cutoff
        :return: Config data
        """
        config_data = {
            'sample_name': sample_name,
            'working_dir': self._working_dir,
            'assembly': {'spades': {}}
        }

        # SPAdes options
        if kmers is not None:
            config_data['assembly']['spades']['kmers'] = kmers
        if cov_cutoff is not None:
            config_data['assembly']['spades']['cov_cutoff'] = cov_cutoff

        # Length filtering
        if min_contig_length is not None:
            # noinspection PyTypeChecker
            config_data['assembly']['min_contig_length'] = min_contig_length

        return config_data

    def __prepare_input_files(self, reads_pe: List[ToolIOFile], reads_se_fwd: List[ToolIOFile],
                              reads_se_rev: List[ToolIOFile]) -> None:
        """
        Prepares the input files for the assembly workflow.
        :return: None
        """
        dir_trimming_out = os.path.dirname(os.path.join(self._working_dir, OUTPUT_READ_TRIMMING_READS_PE))
        if not os.path.isdir(dir_trimming_out):
            os.makedirs(dir_trimming_out)
        SnakemakeUtils.dump_object(reads_pe, os.path.join(self._working_dir, OUTPUT_READ_TRIMMING_READS_PE))
        SnakemakeUtils.dump_object(reads_se_fwd, os.path.join(self._working_dir, OUTPUT_READ_TRIMMING_READS_SE_FWD))
        SnakemakeUtils.dump_object(reads_se_rev, os.path.join(self._working_dir, OUTPUT_READ_TRIMMING_READS_SE_REV))

    def __set_output(self, output_files: Dict[str, str]) -> None:
        """
        Runs the Snakemake workflow.
        :param output_files: Output files dictionary
        :return: None
        """
        log_file_path = os.path.join(self._working_dir, 'camel.log')
        informs = [SnakemakeUtils.load_object(output_files['INFORMS_spades'])]
        if 'INFORMS_seqtk' in output_files:
            informs.append(SnakemakeUtils.load_object(output_files['INFORMS_seqtk']))
        self._output = AssemblyOutput(
            report_section=SnakemakeUtils.load_object(output_files['HTML'])[0].value,
            fasta_contigs=SnakemakeUtils.load_object(output_files['FASTA'])[0],
            informs=informs,
            log_file=log_file_path if os.path.isfile(log_file_path) else None
        )

    @property
    def output(self) -> AssemblyOutput:
        """
        Returns the output of the assembly workflow.
        :return: Assembly output
        """
        return self._output
