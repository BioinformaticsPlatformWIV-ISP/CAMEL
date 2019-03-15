import argparse
import logging
from dataclasses import dataclass
from typing import List, Optional, Dict

import os

from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.workflows.assemblywrapper import AssemblyWrapper
from camel.app.components.workflows.readtrimmingwrapper import ReadTrimmingWrapper
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils


@dataclass
class ReadInput:
    pe: List[ToolIOFile]
    se_fwd: List[ToolIOFile]
    se_rev: List[ToolIOFile]


class MainScriptHelper(object):
    """
    This class is used as helper class for common operation for main scripts.
    For example the read trimming, de-novo assembly that is shared by tools such as 'MLST tool', 'Gene detection tool',
    'ResFinder local' can be done by this class.
    """

    def __init__(self, working_dir: str, sample_name: str):
        """
        Initializes the helper class.
        :param working_dir: Working directory
        """
        self._working_dir = working_dir
        self._sample_name = sample_name
        self._log_files = {}
        self._informs = []

    def trim_reads(self, fastq_reads_raw: List[str], report: Optional[HtmlReport] = None, threads: int = 8,
                   export_fastq: bool = False) -> ReadInput:
        """
        Runs the read trimming workflow.
        :param fastq_reads_raw: Raw FASTQ PE reads
        :param report: If set, the output is added to the given report
        :param threads: Number of threads
        :param export_fastq: If True, FASTQ files are included in the report
        :return: Tuple of trimmed PE reads, trimmed forward reads, trimmed reverse reads
        """
        logging.info("Trimming FASTQ input reads")
        trimming = ReadTrimmingWrapper(os.path.join(self._working_dir, 'trimming'))
        trimming.run_workflow(fastq_reads_raw, threads, export_fastq)
        if report is not None:
            report.add_html_object(trimming.output.report_section)
            trimming.output.report_section.copy_files(report.output_dir)
            report.save()
        if trimming.output.log_file is not None:
            self._log_files['trimming'] = trimming.output.log_file
        self._informs.append(trimming.output.informs_trimmomatic)
        return ReadInput(trimming.output.trimmed_reads_pe, trimming.output.trimmed_reads_se_fwd,
                         trimming.output.trimmed_reads_se_rev)

    def symlink_fastq_pe_input(self, fastq_pe: List[str], fastq_names: List[str], working_dir: str) -> ReadInput:
        """
        Symlinks the FASTQ PE input.
        :return: Assembly input object
        """
        logging.info("Symlinking FASTQ input reads in working directory")
        links = SnakePipelineUtils.symlink_input_files(
            os.path.join(working_dir, 'input'), fastq_pe, fastq_names)
        return ReadInput([ToolIOFile(x) for x in links], [], [])

    def assemble_fastq_reads(self, assembly_input: ReadInput, report: Optional[HtmlReport] = None,
                             kmers: Optional[str] = None, threads: int = 8) -> ToolIOFile:
        """
        Assembles FASTQ reads using SPAdes
        :param assembly_input: Assembly input
        :param report: If set, the output is added to the given report
        :param kmers: Comma separated list of Kmers to use for the assembly
        :param threads: Number of threads to use
        :return: ToolIOFile FASTA object with the assembled contigs
        """
        logging.info("Starting de-novo assembly")
        assembly = AssemblyWrapper(os.path.join(self._working_dir, 'assembly'))
        assembly.run_workflow(
            self._sample_name, assembly_input.pe, assembly_input.se_fwd, assembly_input.se_rev, kmers, threads)
        if report is not None:
            report.add_html_object(assembly.output.report_section)
            assembly.output.report_section.copy_files(report.output_dir)
            report.save()
        if assembly.output.log_file is not None:
            self._log_files['assembly'] = assembly.output.log_file
        self._informs.append(assembly.output.informs)
        return assembly.output.fasta_contigs

    @staticmethod
    def determine_sample_name(args: argparse.Namespace) -> str:
        """
        Determines the sample names based on the given command line arguments.
        :return: Sample name
        """
        if ('sample_name' in args) and (args.sample_name is not None):
            return args.sample_name
        elif ('fasta_name' in args) and (args.fasta_name is not None):
            return os.path.splitext(args.fasta_name)[0]
        elif ('fasta' in args) and (args.fasta is not None):
            return os.path.splitext(os.path.basename(args.fasta))[0]
        elif args.fastq_pe_names is not None:
            try:
                return FastqUtils.get_sample_name(args.fastq_pe_names[0])
            except ValueError as e:
                logging.warning(str(e))
        elif args.fastq_pe is not None:
            try:
                return FastqUtils.get_sample_name(args.fastq_pe[0])
            except ValueError as e:
                logging.warning(str(e))
        logging.warning("Cannot determine sample name from given arguments")
        return 'NA'

    @staticmethod
    def determine_input_files(args: argparse.Namespace) -> str:
        """
        Determines the input files based on the given command line arguments.
        :param args: Command line arguments
        :return: Input files as string
        """
        if (args.fasta is not None) and (args.fasta_name is not None):
            return args.fasta_name
        elif args.fasta is not None:
            return os.path.basename(args.fasta)
        elif (args.fastq_pe is not None) and (args.fastq_pe_names is not None):
            return ', '.join(args.fastq_pe_names)
        elif args.fastq_pe is not None:
            return ', '.join([os.path.basename(f) for f in args.fastq_pe])
        logging.warning("Cannot determine input files from given arguments")
        return 'NA'

    @staticmethod
    def add_common_arguments(argument_parser: argparse.ArgumentParser) -> None:
        """
        Adds the common arguments to the argument parser.
        :param argument_parser: Argument parser
        :return: None
        """
        argument_parser.add_argument('--sample-name', type=str)
        argument_parser.add_argument('--output-dir', required=True, type=str)
        argument_parser.add_argument('--output-html', required=True, type=str)
        argument_parser.add_argument('--working-dir', default=os.path.abspath('.'), type=str)
        argument_parser.add_argument('--threads', default=8, type=int)

    @staticmethod
    def prepare_galaxy_output(output_dir: str, output_html: str) -> None:
        """
        Prepares the Galaxy output files at the start of the script.
        - The output HTML file is removed, so Snakemake can regenerate it
        - The output directory is created if it does not exist yet.
        :param output_dir: Output directory
        :param output_html: Output report path
        :return: None
        """
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        if os.path.isfile(output_html):
            os.remove(output_html)

    @property
    def logs(self) -> Dict[str, str]:
        """
        Returns the log files (key: name, value: log file path).
        :return: Logs
        """
        return self._log_files

    @property
    def informs(self) -> List[Dict[str, str]]:
        """
        Returns the informs.
        :return: List of informs
        """
        return self._informs

    def get_blast_input(self, args: argparse.Namespace, report: Optional[HtmlReport] = None) -> ToolIOFile:
        """
        Returns the input for BLAST based detection methods.
        Takes into accounts:
        - Type of input (FASTA / FASTQ)
        - Read trimming (optional)
        - Kmer setting for de-novo assembly
        :param args: Command line arguments
        :param report: HTML report
        :return: FASTA input for BLAST.
        """
        dir_input = os.path.join(self._working_dir, 'input')

        # FASTA input
        if args.fasta is not None:
            if args.fasta_name is not None:
                link = SnakePipelineUtils.symlink_input_files(dir_input, [args.fasta], [args.fasta_name], True)[0]
                return ToolIOFile(link)
            else:
                return ToolIOFile(args.fasta)

        # FASTQ input
        names = args.fastq_pe_names if args.fastq_pe_names else [os.path.basename(f) for f in args.fastq_pe]
        fq_links = SnakePipelineUtils.symlink_input_files(dir_input, args.fastq_pe, names, True)
        if args.trim_reads:
            assembly_input = self.trim_reads(fq_links, report, args.threads, args.report_include_fastq)
        else:
            assembly_input = ReadInput([ToolIOFile(l) for l in fq_links], [], [])
        return self.assemble_fastq_reads(assembly_input, report, args.kmers, args.threads)

    def get_srst2_input(self, args: argparse.Namespace, report: Optional[HtmlReport] = None) -> List[ToolIOFile]:
        """
        Returns the input for SRST2 (forward + reverse reads).
        :param args: Command line arguments
        :param report: HTML report
        :return: SRST2 input files
        """
        dir_input = os.path.join(self._working_dir, 'input')
        names = args.fastq_pe_names if args.fastq_pe_names else [os.path.basename(f) for f in args.fastq_pe]
        fq_links = SnakePipelineUtils.symlink_input_files(dir_input, args.fastq_pe, names, True)
        if args.trim_reads is False:
            return [ToolIOFile(l) for l in fq_links]
        trimming_out = self.trim_reads(fq_links, report, args.threads, args.report_include_fastq)
        return trimming_out.pe
