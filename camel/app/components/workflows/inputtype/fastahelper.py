import argparse
from pathlib import Path

from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.workflows.inputtype.inputtypehelperbase import InputTypeHelperBase
from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.loggers import logger


class FastaHelper(InputTypeHelperBase):
    """
    Helper class for FASTA input.
    """

    def trim_reads(self, fastq_input: FastqInput, report: HtmlReport, include_fastq: bool, threads: int, **kwargs) -> \
            FastqInput:
        """
        Base function for read-type specific trimming.
        :return: Trimmed FASTQ files
        """
        raise NotImplementedError('Filtering is not implemented for FASTA input')

    def prepare_fasta_input(self, report: HtmlReport, args: argparse.Namespace) -> Path:
        """
        Prepares the FASTA input.
        :param report: HTML report
        :param args: Command-line arguments
        :return: FASTA file
        """
        logger.info("Preparing FASTA input (FASTA input)")
        fasta_file = self.symlink_input_files([Path(args.fasta)], [args.fasta_name])[0]
        return Path(fasta_file)
