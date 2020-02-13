from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Dict, Any, Union

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import assembly_spades


@dataclass
class AssemblyOutput:
    report_section: HtmlReportSection
    fasta_contigs: ToolIOFile
    log_file: Optional[Path]
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
        self._working_dir = Path(working_dir)
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
        if not self._working_dir.exists():
            self._working_dir.mkdir(parents=True)
        self.__prepare_input_files(reads_pe, reads_se_fwd, reads_se_rev)
        config_data = self.__get_config_data(sample_name, kmers, cov_cutoff, min_contig_length)
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._working_dir)
        output_files = {
            'HTML': self._working_dir / assembly_spades.OUTPUT_ASSEMBLY_REPORT,
            'FASTA': self._working_dir / assembly_spades.OUTPUT_ASSEMBLY_FASTA,
            'INFORMS_spades': self._working_dir / assembly_spades.OUTPUT_ASSEMBLY_INFORMS
        }
        if min_contig_length is not None:
            output_files['INFORMS_seqtk'] = self._working_dir / assembly_spades.OUTPUT_ASSEMBLY_FILTERING_INFORMS
        SnakePipelineUtils.run_snakemake(
            assembly_spades.SNAKEFILE_ASSEMBLY_SPADES, config_file, list(output_files.values()), self._working_dir,
            threads)
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
        fq_dict = {'PE': reads_pe, 'SE_FWD': reads_se_fwd, 'SE_REV': reads_se_rev}
        SnakemakeUtils.dump_object(fq_dict, self._working_dir / 'fq_dict.io')

    def __set_output(self, output_files: Dict[str, str]) -> None:
        """
        Runs the Snakemake workflow.
        :param output_files: Output files dictionary
        :return: None
        """
        log_file_path = self._working_dir / 'camel.log'
        informs = [SnakemakeUtils.load_object(output_files['INFORMS_spades'])]
        if 'INFORMS_seqtk' in output_files:
            informs.append(SnakemakeUtils.load_object(output_files['INFORMS_seqtk']))
        self._output = AssemblyOutput(
            report_section=SnakemakeUtils.load_object(output_files['HTML'])[0].value,
            fasta_contigs=SnakemakeUtils.load_object(output_files['FASTA'])[0],
            informs=informs,
            log_file=log_file_path if log_file_path.exists() else None
        )

    @property
    def output(self) -> AssemblyOutput:
        """
        Returns the output of the assembly workflow.
        :return: Assembly output
        """
        return self._output
