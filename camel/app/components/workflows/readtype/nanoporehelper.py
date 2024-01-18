import argparse
from pathlib import Path
from typing import Union, List, Dict

import logging

from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.workflows.readtype.basereadtypehelper import BaseReadTypeHelper
from camel.app.components.workflows.trimmingontwrapper import TrimmingONTWrapper
from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.io.toolio import ToolIO
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.seqtk.seqtkconvert import SeqtkConvert

class NanoporeHelper(BaseReadTypeHelper):
    """
    Helper class for Nanopore reads.
    """

    def __symlink_nanopore_reads(self, fastq_file: Union[Path, None], sample_name: str) -> Path:
        """
        Symlinks the input files to a standardized format based on the sample name.
        :param fastq_file: Input FASTQ file
        :param sample_name: Sample name
        :return: Path to renamed file
        """
        if fastq_file is None:
            raise ValueError("Nanopore data should be SE")
        new_name = f"{sample_name}.fastq{'.gz' if FileSystemHelper.is_gzipped(fastq_file) else ''}"
        return self.symlink_input_files([Path(fastq_file)], [new_name])[0]

    def trim_reads(self, fastq_input: FastqInput, report: HtmlReport, include_fastq: bool, threads: int = 4,
                   **kwargs) -> FastqInput:
        """
        Trims Nanopore reads using filtlong.
        :param fastq_input: FASTQ input
        :param report: HTML report
        :param include_fastq: Boolean to indicate if FASTQ files should be included in the report
        :param threads: Nb. of threads
        :return: FastqInput object with trimmed reads
        """
        logger.info("Trimming reads (Nanopore data)")
        if fastq_input.is_pe or fastq_input.se is None:
            raise ValueError("Nanopore input should be SE")
        # Run workflow
        trimming = TrimmingONTWrapper(self._working_dir / 'trimming')
        trimming.run_workflow(Path(fastq_input.se[0].path), include_fastq, threads)

        # Save output
        report_section_trimming = trimming.output.report_section
        report.add_html_object(report_section_trimming)
        report_section_trimming.copy_files(report.output_dir)
        self._log_files['trimming'] = trimming.output.log_file
        return FastqInput('nanopore', se=trimming.output.trimmed_reads, is_pe=False)

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
        fq_input_se = self.__symlink_nanopore_reads(Path(args.fastq_se), self._sample_name)
        fastq_input = FastqInput(args.read_type, se=[ToolIOFile(Path(fq_input_se))], is_pe=False)
        if args.trim_reads:
            return self.trim_reads(fastq_input, report, args.report_include_fastq, args.threads)
        else:
            return FastqInput(args.read_type, se=[ToolIOFile(Path(fq_input_se))], is_pe=False)

    def prepare_fasta_read_input(self, report: HtmlReport, args: argparse.Namespace) -> dict[
        str, list[Union[ToolIOFile, ToolIOValue, ToolIODirectory, ToolIO]]]:
        """
        Prepares the FASTA input by converting FASTQ to FASTA with Seqtk
        :param report: HTML report
        :param args: Command-line arguments
        :return: FASTA file
        """
        fq_input_se = self.__symlink_nanopore_reads(Path(args.fastq_se), self._sample_name)
        fastq_input = FastqInput(args.read_type, se=[ToolIOFile(Path(fq_input_se))], is_pe=False)
        if args.trim_reads:
            convert_input = self.trim_reads(fastq_input, report, args.report_include_fastq, args.threads)
        else:
            convert_input = fastq_input
        camel = Camel.get_instance()
        convert = SeqtkConvert(camel)
        convert.add_input_files({'FASTQ': convert_input.se})
        convert.run(self._working_dir)
        return Path(convert._output_string)