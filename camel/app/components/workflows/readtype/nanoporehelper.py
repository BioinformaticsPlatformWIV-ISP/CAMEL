import argparse
from pathlib import Path

from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.workflows.readtype.basereadtypehelper import BaseReadTypeHelper
from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.io.tooliofile import ToolIOFile


class NanoporeHelper(BaseReadTypeHelper):
    """
    Helper class for Nanopore reads.
    """

    def trim_reads(self, fastq_input: FastqInput, report: HtmlReport, include_fastq: bool, threads: int) -> FastqInput:
        """
        Trims Illumina reads using Trimmomatic.
        :param fastq_input: FASTQ input
        :param report: HTML report
        :param include_fastq: Boolean to indicate if FASTQ files should be included in the report
        :param threads: Nb. of threads
        :return: FastqInput object with trimmed reads
        """
        raise NotImplementedError("Trimming for Nanopore data is currently not supported")

    def prepare_fasta_input(self, report: HtmlReport, args: argparse.Namespace) -> Path:
        """
        Prepares the FASTA input.
        :param report: HTML report
        :param args: Command-line arguments
        :return: FASTA file
        """
        raise NotImplementedError("Assembly for Nanopore data is currently not supported.")

    def prepare_fastq_input(self, report: HtmlReport, args: argparse.Namespace) -> FastqInput:
        """
        Prepares the FASTQ input.
        :param report: HTML report
        :param args: Command-line arguments
        :return: FASTQ input
        """
        fq_input_se = self.symlink_input_files([Path(args.fastq_se)], [args.fastq_se_name])[0]
        fastq_input = FastqInput(args.read_type, se=[ToolIOFile(fq_input_se)], is_pe=False)
        if args.trim_reads:
            return self.trim_reads(fastq_input, report, args.report_include_fastq, args.threads)
        else:
            return FastqInput(args.read_type, se=[ToolIOFile(fq_input_se)], is_pe=False)
