from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Any, Dict

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.sequencetyping.sequencetypingutils import SequenceTypingUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import sequence_typing, assembly_spades


@dataclass
class SequenceTypingInput:
    """
    This class is used to construct the input for the sequence typing workflow.
    """
    sample_name: str
    db_path: str
    fasta: Optional[ToolIOFile] = None
    fastq_pe: Optional[List[ToolIOFile]] = None
    db_key: str = 'mlst'


@dataclass
class SequenceTypingOutput:
    """
    This class contains the output of the sequence typing workflow.
    """
    report_section: HtmlReportSection
    informs: List[Dict[str, Any]]
    log_file: Optional[str] = None


class SequenceTypingWrapper(object):
    """
    This class is used as a wrapper class around the sequence typing Snakemake workflow.
    """

    def __init__(self, working_dir: str) -> None:
        """
        Initializes the read sequence typing helper.
        :param working_dir: Working directory
        """
        self._working_dir = Path(working_dir)
        self._output = None

    def run_workflow_blast(self, workflow_input: SequenceTypingInput, blast_task: Optional[str] = None,
                           threads: Optional[int] = 8) -> None:
        """
        Runs the BLAST based sequence typing workflow.
        :param workflow_input: Workflow input
        :param blast_task: blast task
        :param threads: Number of threads to use
        :return: None
        """
        if not self._working_dir.exists():
            self._working_dir.mkdir(parents=True)
        if workflow_input.fasta is None:
            raise ValueError("FASTA input is required to run the blast analysis")
        self.__create_blast_input(workflow_input.fasta.path)
        options = {'blastn_task': blast_task} if blast_task is not None else None
        config_path = self.__create_config_file(
            workflow_input.sample_name, 'blast', workflow_input.db_key, workflow_input.db_path, options)
        self.__run_snakefile(config_path, workflow_input.db_key, threads)

    def run_workflow_srst2(self, workflow_input: SequenceTypingInput, srst2_options: Optional[Dict[str, Any]] = None,
                           threads: Optional[int] = 8) -> None:
        """
        Runs the SRST2 based sequence typing workflow.
        :param workflow_input: Input for the workflow
        :param srst2_options: Additional options for SRST2
        :param threads: Number of threads to use
        :return: None
        """
        if not self._working_dir.exists():
            self._working_dir.mkdir()
        if workflow_input.fastq_pe is None:
            raise ValueError('FASTQ input is required to run the SRST2 analysis')
        scheme_metadata = SequenceTypingUtils.parse_scheme_metadata(workflow_input.db_path)
        if len([locus for locus in scheme_metadata['loci'] if locus['type'] == 'peptide']) > 0:
            if workflow_input.fasta is None:
                raise ValueError('FASTA input is required when there are protein loci in the scheme')
            else:
                self.__create_blast_input(workflow_input.fasta.path)
        self.__create_srst2_input([f.path for f in workflow_input.fastq_pe])
        additional_options = srst2_options if srst2_options else None
        config_path = self.__create_config_file(
            workflow_input.sample_name, 'srst2', workflow_input.db_key, workflow_input.db_path, additional_options)
        self.__run_snakefile(config_path, workflow_input.db_key, threads)

    def run_workflow_kma(self, workflow_input: SequenceTypingInput, threads: Optional[int] = 8) -> None:
        """
        Runs the KMA based sequence typing workflow.
        :param workflow_input: Input for the workflow
        :param threads: Number of threads to use
        :return: None
        """
        if not self._working_dir.exists():
            self._working_dir.mkdir()
        if workflow_input.fastq_pe is None:
            raise ValueError('FASTQ input is required to run the KMA analysis')
        scheme_metadata = SequenceTypingUtils.parse_scheme_metadata(workflow_input.db_path)
        if len([locus for locus in scheme_metadata['loci'] if locus['type'] == 'peptide']) > 0:
            if workflow_input.fasta is None:
                raise ValueError('FASTA input is required when there are protein loci in the scheme')
            else:
                self.__create_blast_input(workflow_input.fasta.path)
        self.__create_srst2_input([f.path for f in workflow_input.fastq_pe])
        config_path = self.__create_config_file(
            workflow_input.sample_name, 'kma', workflow_input.db_key, workflow_input.db_path, {})
        self.__run_snakefile(config_path, workflow_input.db_key, threads)

    def __create_config_file(self, sample_name: str, detection_method: str, db_key: str, db_path: str,
                             additional_options: Dict[str, Any]) -> str:
        """
        Creates the configuration file for the sequence typing workflow.
        :param sample_name: Sample name
        :param detection_method: Detection method ('blast', 'srst2)
        :param db_key: Database key (e.g. 'mlst')
        :param db_path: Database path
        :param additional_options: Additional options to add to the configuration file.
        :return: Config file path
        """
        data = {
            'working_dir': str(self._working_dir),
            'sample_name': sample_name,
            'detection_method': detection_method,
            'sequence_typing': {
                db_key: {
                    'path': db_path,
                    **additional_options
                }
            }
        }
        import pprint
        pprint.pprint(data)
        return SnakePipelineUtils.generate_config_file(data, self._working_dir)

    def __create_blast_input(self, fasta_path: str) -> None:
        """
        Creates the input for the workflow.
        :param fasta_path: FASTA file path
        :return: None
        """
        fasta_path_output = (self._working_dir / assembly_spades.OUTPUT_ASSEMBLY_FASTA)
        if not fasta_path_output.parent.exists():
            fasta_path_output.parent.mkdir(parents=True)
        SnakemakeUtils.dump_object([ToolIOFile(fasta_path)], str(fasta_path_output))

    def __create_srst2_input(self, fastq_pe: List[str]) -> None:
        """
        Creates the input for the SRST2 workflow.
        :param fastq_pe: Paired-end FASTQ files
        :return: None
        """
        fq_dict = {'PE': [ToolIOFile(x) for x in fastq_pe]}
        SnakemakeUtils.dump_object(fq_dict, str(self._working_dir / 'fq_dict.io'))

    def __run_snakefile(self, config_path: str, db_key: str, threads: int) -> None:
        """
        Runs the Snakefile.
        :param config_path: Path to the configuration file
        :param db_key: Database key
        :param threads: Number of threads to use
        :return: None
        """
        output_files = {
            'report': self._working_dir / str(sequence_typing.OUTPUT_TYPING_REPORT).format(scheme=db_key),
            'informs': self._working_dir / str(sequence_typing.OUTPUT_TYPING_INFORMS).format(scheme=db_key)
        }
        # Run Snakemake
        SnakePipelineUtils.run_snakemake(
            sequence_typing.SNAKEFILE_SEQUENCE_TYPING, config_path, list(output_files.values()), self._working_dir,
            threads)

        # Collect output
        log_file_path = self._working_dir / 'camel.log'
        self._output = SequenceTypingOutput(
            report_section=SnakemakeUtils.load_object(str(output_files['report']))[0].value,
            informs=SnakemakeUtils.load_object(str(output_files['informs'])),
            log_file=log_file_path if log_file_path.exists() else None
        )

    @property
    def output(self) -> SequenceTypingOutput:
        """
        Returns the output generated by the gene detection workflow.
        :return: Output
        """
        return self._output
