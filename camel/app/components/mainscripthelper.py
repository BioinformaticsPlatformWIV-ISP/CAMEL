import argparse
import logging
from typing import List, Optional, Tuple

import os

from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.workflows.assemblywrapper import AssemblyWrapper
from camel.app.components.workflows.readtrimmingwrapper import ReadTrimmingWrapper
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils


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

    def trim_reads(self, fastq_reads_raw: List[str], report: Optional[HtmlReport]=None) -> \
            Tuple[List[ToolIOFile], List[ToolIOFile], List[ToolIOFile]]:
        """
        Runs the read trimming workflow.
        :param fastq_reads_raw: Raw FASTQ PE reads
        :param report: If set, the output is added to the given report
        :return: Tuple of trimmed PE reads, trimmed forward reads, trimmed reverse reads
        """
        trimming = ReadTrimmingWrapper(os.path.join(self._working_dir, 'trimming'))
        trimming.run_workflow(fastq_reads_raw)
        if report is not None:
            report.add_html_object(trimming.output.report_section)
            trimming.output.report_section.copy_files(report.output_dir)
            report.save()
        return trimming.output.trimmed_reads_pe, trimming.output.trimmed_reads_se_fwd, \
            trimming.output.trimmed_reads_se_rev

    def assemble_fastq_reads(self, fastq_pe: List[str], fastq_names: List[str], perform_trimming: bool,
                             report: Optional[HtmlReport]=None, kmers: Optional[str]=None) -> ToolIOFile:
        """
        Assembles FASTQ reads using SPAdes
        :param fastq_pe: FASTQ PE read paths
        :param fastq_names: FASTQ file names
        :param perform_trimming: If True, reads are trimmed before the assembly
        :param report: If set, the output is added to the given report
        :param kmers: Comma separated list of Kmers to use for the assembly
        :return: ToolIOFile FASTA object with the assembled contigs
        """
        if perform_trimming is True:
            assembly_input = self.trim_reads(fastq_pe, report)
        else:
            links = SnakePipelineUtils.symlink_input_files(
                os.path.join(self._working_dir, 'input'), fastq_pe, fastq_names)
            assembly_input = ([ToolIOFile(x) for x in links], [], [])
        assembly = AssemblyWrapper(os.path.join(self._working_dir, 'assembly'))
        assembly.run_workflow(self._sample_name, assembly_input[0], assembly_input[1], assembly_input[2], kmers)
        if report is not None:
            report.add_html_object(assembly.output.report_section)
            assembly.output.report_section.copy_files(report.output_dir)
            report.save()
        return assembly.output.fasta_contigs

    @staticmethod
    def determine_sample_name(args: argparse.Namespace) -> str:
        """
        Determines the sample names based on the given command line arguments.
        :return: Sample name
        """
        if args.sample_name is not None:
            return args.sample_name
        elif args.fasta_name is not None:
            return os.path.splitext(args.fasta_name)[0]
        elif args.fasta is not None:
            return os.path.splitext(os.path.basename(args.fasta_path))[0]
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
