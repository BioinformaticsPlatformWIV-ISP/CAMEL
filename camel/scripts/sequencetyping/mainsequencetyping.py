#!/usr/bin/env python
import argparse
import datetime
import json
import logging
from typing import Any, Dict, List

import os

from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.mainscripthelper import MainScriptHelper
from camel.app.components.workflows.sequencetypingwrapper import SequenceTypingWrapper
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources import CSS_STYLE
from camel.resources.javascript import JQUERY_SRC


class MainSequenceTyping(object):
    """
    Class to run sequence typing tool, it supports both BLAST+ and SRST2 as detection methods for alleles.
    """

    def __init__(self, args: argparse.Namespace = None):
        """
        Initializes the main script.
        """
        self._args = args if args is not None else MainSequenceTyping._parse_arguments()
        self._sample_name = MainScriptHelper.determine_sample_name(self._args)
        self._helper = MainScriptHelper(self._args.working_dir, self._sample_name)
        self._report = None

    @staticmethod
    def _parse_arguments() -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        MainScriptHelper.add_common_arguments(argument_parser)
        argument_parser.add_argument('--fasta', help="Input FASTA file", type=str)
        argument_parser.add_argument('--fasta-name', help="Input FASTA file name", type=str)
        argument_parser.add_argument('--fastq-pe', help="Input PE FASTQ files", nargs=2)
        argument_parser.add_argument('--fastq-pe-names', help="Input PE FASTQ file names", nargs=2)
        argument_parser.add_argument('--trim-reads', help="Perform read trimming", action='store_true')
        argument_parser.add_argument('--kmers', help="Kmers to use for assembly", type=str)
        argument_parser.add_argument('--scheme-dir', required=True, type=str)
        argument_parser.add_argument('--detection-method', type=str, choices=['blast', 'srst2'], default='blast')
        argument_parser.add_argument('--report-include-fastq', action='store_true')
        return argument_parser.parse_args()

    def run(self):
        """
        Runs the workflow.
        :return: None
        """
        self.__init_report()
        self.__add_analysis_info_section()
        db_data = self.__get_db_metadata(self._args.scheme_dir)
        if self._args.detection_method == 'blast':
            fasta_file = self.__get_fasta_file_blast()
            self.__run_sequence_typing_blast(fasta_file.path, db_data['name'], self._args.scheme_dir)
        elif self._args.detection_method == 'srst2':
            input_pe = self.__get_paired_reads_srst2()
            self.__run_sequence_typing_srst2(input_pe, db_data['name'], self._args.scheme_dir)

    def __init_report(self) -> None:
        """
        Initializes the HTML report
        :return: None
        """
        self._report = HtmlReport(self._args.output_html, self._args.output_dir, [JQUERY_SRC])
        if not os.path.isdir(self._args.output_dir):
            os.makedirs(self._args.output_dir)
        self._report.initialize('Sequence typing report', CSS_STYLE)
        self._report.add_pipeline_header(f'Sequence typing ({self._args.detection_method})')
        self._report.save()

    def __add_analysis_info_section(self) -> None:
        """
        Adds the report section with the analysis info
        :return: None
        """
        section = HtmlReportSection('Analysis info')
        section.add_table([
            ['Analysis date:', datetime.datetime.now().strftime(SnakePipelineUtils.DATE_FORMAT)],
            ['Input file(s):', self._helper.determine_input_files(self._args)],
        ], table_attributes=[('class', 'information')])
        self._report.add_html_object(section)
        self._report.save()

    def __get_fasta_file_blast(self) -> ToolIOFile:
        """
        Returns the input FASTA file
        :return: FASTA file
        """
        if self._args.fasta is not None:
            return ToolIOFile(self._args.fasta)
        else:
            if self._args.trim_reads:
                assembly_input = self._helper.trim_reads(
                    self._args.fastq_pe, self._report, self._args.threads, self._args.report_include_fastq)
            else:
                assembly_input = self._helper.symlink_fastq_pe_input(
                    self._args.fastq_pe, self._args.fastq_pe_names, self._args.working_dir)
            return self._helper.assemble_fastq_reads(assembly_input, self._report, self._args.kmers, self._args.threads)

    def __get_paired_reads_srst2(self) -> List[str]:
        """
        Returns the FASTQ PE input for SRST2.
        :return: FASTQ PE files
        """
        if self._args.fastq_pe is None:
            raise ValueError("FASTQ PE input needs to be provided")
        if self._args.trim_reads is True:
            trimming_output = self._helper.trim_reads(
                self._args.fastq_pe, self._report, self._args.threads, self._args.report_include_fastq)
            return [x.path for x in trimming_output.pe]
        else:
            gzipped = self._args.fastq_pe[0].endswith('.gz')
            return SnakePipelineUtils.symlink_input_files(
                os.path.join(self._args.working_dir, 'input'), self._args.fastq_pe,
                [f"{self._sample_name}_{x}.fastq{'.gz' if gzipped else ''}" for x in (1, 2)])

    def __get_db_metadata(self, directory: str) -> Dict[str, Any]:
        """
        Returns the database metadata.
        :param directory: Database directory
        :return: Metadata
        """
        with open(os.path.join(directory, 'scheme_metadata.txt')) as handle:
            return json.load(handle)

    def __run_sequence_typing_blast(self, fasta_file: str, db_key: str, db_path: str) -> None:
        """
        Runs the sequence typing workflow using BLAST.
        :param fasta_file: Input FASTA file
        :param db_key: Database key
        :param db_path: Database directory path
        :return: None
        """
        wrapper = SequenceTypingWrapper(self._args.working_dir)
        workflow_input = SequenceTypingWrapper.SequenceTypingInput(
            sample_name=self._sample_name, fasta=ToolIOFile(fasta_file), db_path=db_path, db_key=db_key)
        wrapper.run_workflow_blast(workflow_input, self._args.threads)
        self._report.add_html_object(wrapper.output.report_section)
        wrapper.output.report_section.copy_files(self._report.output_dir)
        self._report.save()

    def __run_sequence_typing_srst2(self, fastq_pe: List[str], db_key: str, db_path: str) -> None:
        """
        Runs the sequence typing workflow using SRST2.
        :param fastq_pe: In put FASTQ PE files
        :param db_key: Database key
        :param db_path: Database path
        :return: None
        """
        wrapper = SequenceTypingWrapper(self._args.working_dir)
        workflow_input = SequenceTypingWrapper.SequenceTypingInput(
            fasta=ToolIOFile(self._args.fasta) if self._args.fasta else None,
            sample_name=self._sample_name, fastq_pe=[ToolIOFile(x) for x in fastq_pe],
            db_key=db_key, db_path=db_path)
        wrapper.run_workflow_srst2(workflow_input, self._args.threads)
        self._report.add_html_object(wrapper.output.report_section)
        wrapper.output.report_section.copy_files(self._report.output_dir)
        self._report.save()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    sequence_typing = MainSequenceTyping()
    sequence_typing.run()
