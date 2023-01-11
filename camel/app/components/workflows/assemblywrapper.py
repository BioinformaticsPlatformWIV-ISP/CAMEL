from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Dict, Any, Union

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import assembly_spades, assembly_canu


@dataclass
class AssemblyOutput:
    report_section: HtmlReportSection
    tsv_summary: Path
    fasta_contigs: Path
    log_file: Optional[Path]
    informs: List[Dict[str, Any]]
    qc_stats: Optional[Dict[str, Any]] = None


class AssemblyWrapper(object):
    """
    This class is used as a wrapper class around the assembly Snakemake workflow.
    """

    def __init__(self, working_dir: Path, read_type: str = 'illumina') -> None:
        """
        Initializes the read trimming helper.
        :param working_dir: Working directory
        :param read_type: Read type
        :return: None
        """
        self._working_dir = working_dir
        self._working_dir.mkdir(parents=True, exist_ok=True)
        self._output = None
        self._read_type = read_type

    @property
    def assembler_key(self) -> str:
        """
        Returns the key for the assembler.
        :return: Assembler key
        """
        if self._read_type in ('illumina', 'iontorrent'):
            return 'spades'
        elif self._read_type == 'nanopore':
            return 'canu'
        raise ValueError(f'Invalid read type: {self._read_type}')

    def run(self, name: str, fastq_in: FastqInput, min_ctg_len: Optional[int] = None,
            assembler_opts: Optional[Dict] = None, calc_qc_stats: bool = False, threads: int = 8) -> None:
        """
        Runs the assembly workflow for paired-end input.
        :param name: Dataset name
        :param fastq_in: FASTQ input
        :param min_ctg_len: Minimum contig length
        :param assembler_opts: Additional options for the assembler
        :param calc_qc_stats: If True, determines the QC stats
        :param threads: Number of threads
        :return: None
        """
        if fastq_in.is_pe:
            fq_dict = {'PE': fastq_in.pe}
            if fastq_in.se_fwd is not None:
                fq_dict['SE_FWD'] = fastq_in.se_fwd
            if fastq_in.se_rev is not None:
                fq_dict['SE_REV'] = fastq_in.se_rev
        else:
            fq_dict = {'SE': fastq_in.se}
        SnakemakeUtils.dump_object(fq_dict, self._working_dir / 'fq_dict.io')
        self.__run_workflow(name, min_ctg_len, assembler_opts, calc_qc_stats, threads)

    def __run_workflow(
            self, name: str, min_ctg_len: Optional[int] = None, assembler_opts: Optional[Dict] = None,
            calc_qc_stats: bool = False, threads: int = 8) -> None:
        """
        Runs the underlying workflow.
        :return: None
        """
        # Creat working directory
        if not self._working_dir.exists():
            self._working_dir.mkdir(parents=True)

        # Create config file
        config_data = self.__get_config_data(name, min_ctg_len, assembler_opts)
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._working_dir)

        # Collect output files
        output_files = self.__get_output_files_dict(min_ctg_len, calc_qc_stats)
        path_workflow = assembly_spades.SNAKEFILE_ASSEMBLY_SPADES if self._read_type != 'nanopore' else \
            assembly_canu.SNAKEFILE_ASSEMBLY_CANU
        SnakePipelineUtils.run_snakemake(
            path_workflow, config_file, list(output_files.values()), Path(self._working_dir), threads)
        self.__set_output(output_files)

    def __get_output_files_dict(self, min_ctg_len: Union[int, None], calc_qc_stats: bool) -> Dict[str, Path]:
        """
        Returns the dictionary with output files.
        :param min_ctg_len: Minimum contig length
        :param calc_qc_stats: If True, QC stats are calculated
        :return: Dictionary with output files.
        """
        # SPAdes
        if self._read_type in ('illumina', 'iontorrent'):
            output_files = {
                'HTML': self._working_dir / assembly_spades.OUTPUT_ASSEMBLY_REPORT,
                'TSV': self._working_dir / assembly_spades.OUTPUT_ASSEMBLY_SUMMARY,
                'FASTA': self._working_dir / assembly_spades.OUTPUT_ASSEMBLY_FASTA,
                'INFORMS_spades': self._working_dir / assembly_spades.OUTPUT_ASSEMBLY_INFORMS
            }
            if min_ctg_len is not None:
                output_files['INFORMS_seqtk'] = self._working_dir / assembly_spades.OUTPUT_ASSEMBLY_FILTERING_INFORMS
            if calc_qc_stats is True:
                output_files['INFORMS_bowtie2'] = self._working_dir / assembly_spades.OUTPUT_ASSEMBLY_MAPPING_INFORMS
                output_files['INFORMS_samtools'] = self._working_dir / assembly_spades.OUTPUT_ASSEMBLY_DEPTH_INFORMS
        # CANU
        else:
            output_files = {
                'HTML': self._working_dir / assembly_canu.OUTPUT_ASSEMBLY_REPORT,
                'TSV': self._working_dir / assembly_canu.OUTPUT_ASSEMBLY_SUMMARY,
                'FASTA': self._working_dir / assembly_canu.OUTPUT_ASSEMBLY_FASTA,
                'INFORMS_canu': self._working_dir / assembly_canu.OUTPUT_ASSEMBLY_INFORMS
            }
            if min_ctg_len is not None:
                output_files['INFORMS_seqtk'] = self._working_dir / assembly_canu.OUTPUT_ASSEMBLY_FILTERING_INFORMS
            if calc_qc_stats is True:
                raise NotImplementedError('Statistics for ONT assembly are currently not implemented')
        return output_files

    def __get_config_data(self, name: str, min_ctg_len: Union[int, None], assembler_opts: Optional[Dict] = None) -> \
            Dict[str, Any]:
        """
        Builds the configuration file to run the assembly workflow.
        :param name: Dataset name
        :param min_ctg_len: Minimal contig length
        :param assembler_opts: Additional options for the assembler
        :return: Config data
        """
        config_data = {
            'sample_name': name,
            'working_dir': str(self._working_dir),
            'assembly': {self.assembler_key: {}},
            'read_type': self._read_type
        }

        # Assembler options
        config_data['assembly'][self.assembler_key] = assembler_opts if assembler_opts is not None else {}
        if self._read_type == 'iontorrent':
            config_data['assembly']['spades']['iontorrent'] = None

        # Length filtering
        if min_ctg_len is not None:
            # noinspection PyTypeChecker
            config_data['assembly']['min_contig_length'] = min_ctg_len
        return config_data

    def __set_output(self, output_files: Dict[str, Path]) -> None:
        """
        Runs the Snakemake workflow.
        :param output_files: Output files dictionary
        :return: None
        """
        log_file_path = self._working_dir / 'camel.log'
        informs = [SnakemakeUtils.load_object(output_files[f'INFORMS_{self.assembler_key}'])]
        if 'INFORMS_seqtk' in output_files:
            informs.append(SnakemakeUtils.load_object(output_files['INFORMS_seqtk']))
        if all(key in output_files for key in ('INFORMS_bowtie2', 'INFORMS_samtools')):
            qc_stats = {
                'depth': SnakemakeUtils.load_object(output_files['INFORMS_samtools']),
                'mapping': SnakemakeUtils.load_object(output_files['INFORMS_bowtie2']),
            }
        else:
            qc_stats = None
        self._output = AssemblyOutput(
            report_section=SnakemakeUtils.load_object(output_files['HTML'])[0].value,
            tsv_summary=output_files['TSV'],
            fasta_contigs=SnakemakeUtils.load_object(output_files['FASTA'])[0].path,
            informs=informs,
            log_file=log_file_path if log_file_path.exists() else None,
            qc_stats=qc_stats
        )

    @property
    def output(self) -> AssemblyOutput:
        """
        Returns the output of the assembly workflow.
        :return: Assembly output
        """
        return self._output
