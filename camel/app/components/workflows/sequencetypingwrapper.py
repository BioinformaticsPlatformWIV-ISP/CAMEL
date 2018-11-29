from dataclasses import dataclass
from typing import Optional, List

import os
import yaml

from camel.app.command.command import Command
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.sequencetyping.sequencetypingutils import SequenceTypingUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.resources.snakefile import SNAKEFILE_SEQUENCE_TYPING
from camel.resources.snakefile.assembly_spades import OUTPUT_ASSEMBLY_FASTA
from camel.resources.snakefile.read_trimming import OUTPUT_READ_TRIMMING_READS_PE
from camel.resources.snakefile.sequence_typing import OUTPUT_TYPING_REPORT


class SequenceTypingWrapper(object):
    """
    This class is used as a wrapper class around the sequence typing Snakemake workflow.
    """

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

    def __init__(self, working_dir: str) -> None:
        """
        Initializes the read trimming helper.
        :param working_dir: Working directory
        """
        self._working_dir = working_dir
        self._output = None

    def run_workflow_blast(self, workflow_input: SequenceTypingInput, threads: Optional[int]=8) -> None:
        """
        Runs the BLAST based sequence typing workflow.
        :param workflow_input: Workflow input
        :param threads: Number of threads to use
        :return: None
        """
        if not os.path.isdir(self._working_dir):
            os.makedirs(self._working_dir)
        if workflow_input.fasta is None:
            raise ValueError("FASTA input is required to run the blast analysis")
        self.__create_blast_input(workflow_input.fasta.path)
        config_path = self.__create_config_file(
            workflow_input.sample_name, 'blast', workflow_input.db_key, workflow_input.db_path)
        self.__run_snakefile(config_path, workflow_input.db_key, threads)

    def run_workflow_srst2(self, workflow_input: SequenceTypingInput, threads: Optional[int]=8) -> None:
        """
        RUns the SRST2 based sequence typing workflow.
        :param workflow_input: Input for the workflow
        :param threads: Number of threads to use
        :return: None
        """
        if not os.path.isdir(self._working_dir):
            os.makedirs(self._working_dir)
        if workflow_input.fastq_pe is None:
            raise ValueError('FASTQ input is required to run the SRST2 analysis')
        if len(SequenceTypingUtils.get_loci(workflow_input.db_path)['peptide']) > 0 and workflow_input.fasta is None:
            raise ValueError('FASTA input is required when there are protein loci in the scheme')
        else:
            self.__create_blast_input(workflow_input.fasta.path)
        self.__create_srst2_input([f.path for f in workflow_input.fastq_pe])
        config_path = self.__create_config_file(
            workflow_input.sample_name, 'srst2', workflow_input.db_key, workflow_input.db_path)
        self.__run_snakefile(config_path, workflow_input.db_key, threads)

    def __create_config_file(self, sample_name: str, detection_method: str, db_key: str, db_path: str) -> str:
        """
        Creates the configuration file for the sequence typing workflow.
        :param sample_name: Sample name
        :param detection_method: Detection method ('blast', 'srst2)
        :param db_key: Database key (e.g. 'mlst')
        :param db_path: Database path
        :return: Config file path
        """
        config_path = os.path.join(self._working_dir, 'config_typing.yml')
        with open(config_path, 'w') as handle:
            yaml.dump({
                'working_dir': self._working_dir,
                'sample_name': sample_name,
                'detection_method': detection_method,
                'sequence_typing': {db_key: db_path}
            }, handle)
        return config_path

    def __create_blast_input(self, fasta_path: str) -> None:
        """
        Creates the input for the workflow.
        :param fasta_path: FASTA file path
        :return: None
        """
        target_dir = os.path.dirname(os.path.join(self._working_dir, OUTPUT_ASSEMBLY_FASTA))
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)
        SnakemakeUtils.dump_object([ToolIOFile(fasta_path)], os.path.join(
            self._working_dir, OUTPUT_ASSEMBLY_FASTA))

    def __create_srst2_input(self, fastq_pe: List[str]) -> None:
        """
        Creates the input for the SRST2 workflow.
        :param fastq_pe: Paired-end FASTQ files
        :return: None
        """
        target_dir = os.path.dirname(os.path.join(self._working_dir, OUTPUT_READ_TRIMMING_READS_PE))
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)
        SnakemakeUtils.dump_object([ToolIOFile(x) for x in fastq_pe], os.path.join(
            self._working_dir, OUTPUT_READ_TRIMMING_READS_PE))

    def __run_snakefile(self, config_path: str, db_key: str, threads: int) -> None:
        """
        Runs the Snakefile.
        :param config_path: Path to the configuration file
        :param db_key: Database key
        :param threads: Number of threads to use
        :return: None
        """
        output_path = os.path.join(self._working_dir, OUTPUT_TYPING_REPORT.format(scheme=db_key))
        command = Command('snakemake --snakefile {} --configfile {} {} --cores {}'.format(
            SNAKEFILE_SEQUENCE_TYPING, config_path, output_path, threads
        ))
        command.run_command(self._working_dir)
        self._output = SequenceTypingWrapper.SequenceTypingOutput(
            report_section=SnakemakeUtils.load_object(output_path)[0].value)

    @property
    def output(self) -> SequenceTypingOutput:
        """
        Returns the output generated by the gene detection workflow.
        :return: Output
        """
        return self._output
