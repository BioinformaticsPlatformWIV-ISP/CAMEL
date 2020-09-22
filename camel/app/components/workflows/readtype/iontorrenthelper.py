import logging
from pathlib import Path

import argparse

from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.workflows.readtype.basereadtypehelper import BaseReadTypeHelper
from camel.app.components.workflows.trimmingiontorrentwrapper import TrimmingIonTorrentWrapper
from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.io.tooliofile import ToolIOFile


class IonTorrentHelper(BaseReadTypeHelper):
    """
    Helper class for IonTorrent reads.
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
        logging.info("Trimming reads (IonTorrent data)")
        if fastq_input.is_pe or fastq_input.se is None:
            raise ValueError("IonTorrent input should be SE")
        # Run workflow
        trimming = TrimmingIonTorrentWrapper(self._working_dir / 'trimming')
        if fastq_input.is_pe:
            raise ValueError("PE data not allowed for IonTorrent")
        trimming.run_workflow(Path(fastq_input.se[0].path), include_fastq, threads)

        # Save output
        report_section_trimming = trimming.output.report_section
        report.add_html_object(report_section_trimming)
        report_section_trimming.copy_files(report.output_dir)
        self._informs.extend(trimming.output.informs_seqtk)
        self._log_files['trimming'] = trimming.output.log_file
        return FastqInput('iontorrent', se=trimming.output.trimmed_reads, is_pe=False)

    def prepare_fasta_input(self, report: HtmlReport, args: argparse.Namespace) -> Path:
        """
        Prepares the FASTA input.
        :param report: HTML report
        :param args: Command-line arguments
        :return: FASTA file
        """
        logging.info("Preparing FASTA input (IonTorrent data)")
        if args.fastq_se is None:
            raise ValueError("IonTorrent data should be SE")
        fq_names = [args.fastq_se_name] if args.fastq_se_name is not None else None
        fq_input_se = self.symlink_input_files([Path(args.fastq_se)], fq_names)[0]
        fastq_input = FastqInput(args.read_type, se=[ToolIOFile(fq_input_se)], is_pe=False)
        if args.trim_reads:
            assembly_input = self.trim_reads(fastq_input, report, args.report_include_fastq, args.threads)
        else:
            assembly_input = fastq_input
        return self.assemble_fastq_reads(assembly_input, args, report)

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
