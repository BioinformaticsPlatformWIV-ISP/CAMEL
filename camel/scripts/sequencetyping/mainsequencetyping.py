#!/usr/bin/env python
import argparse
import datetime
import json
import logging
import shutil
from typing import Any, Dict, List

import os

from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.mainscripthelper import MainScriptHelper
from camel.app.components.workflows.sequencetypingwrapper import SequenceTypingWrapper, SequenceTypingInput, \
    SequenceTypingOutput
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
        argument_parser.add_argument('--srst2-max-unaligned-overlap', type=int, default=100)
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
            fasta_file = self._helper.get_blast_input(self._args, self._report)
            output = self.__run_sequence_typing_blast(fasta_file, db_data['name'], self._args.scheme_dir)
        elif self._args.detection_method == 'srst2':
            input_pe = self._helper.get_srst2_input(self._args, self._report)
            output = self.__run_sequence_typing_srst2(input_pe, db_data['name'], self._args.scheme_dir)
        else:
            raise ValueError(f"Invalid detection method: {self._args.detection_method}")
        self.__export_output(output)

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

    def __get_db_metadata(self, directory: str) -> Dict[str, Any]:
        """
        Returns the database metadata.
        :param directory: Database directory
        :return: Metadata
        """
        with open(os.path.join(directory, 'scheme_metadata.txt')) as handle:
            return json.load(handle)

    def __run_sequence_typing_blast(self, fasta_file: ToolIOFile, db_key: str, db_path: str) -> SequenceTypingOutput:
        """
        Runs the sequence typing workflow using BLAST.
        :param fasta_file: Input FASTA file
        :param db_key: Database key
        :param db_path: Database directory path
        :return: None
        """
        wrapper = SequenceTypingWrapper(self._args.working_dir)
        workflow_input = SequenceTypingInput(
            sample_name=self._sample_name, fasta=ToolIOFile(fasta_file.path), db_path=db_path, db_key=db_key)
        wrapper.run_workflow_blast(workflow_input, self._args.threads)
        return wrapper.output

    def __run_sequence_typing_srst2(self, fastq_pe: List[ToolIOFile], db_key: str, db_path: str) -> \
            SequenceTypingOutput:
        """
        Runs the sequence typing workflow using SRST2.
        :param fastq_pe: Input FASTQ PE files
        :param db_key: Database key
        :param db_path: Database path
        :return: None
        """
        wrapper = SequenceTypingWrapper(self._args.working_dir)
        workflow_input = SequenceTypingInput(
            fasta=ToolIOFile(self._args.fasta) if self._args.fasta else None,
            sample_name=self._sample_name, fastq_pe=fastq_pe, db_key=db_key, db_path=db_path)
        srst2_options = {'max_unaligned_overlap': self._args.srst2_max_unaligned_overlap}
        wrapper.run_workflow_srst2(workflow_input, srst2_options, self._args.threads)
        return wrapper.output

    def __export_output(self, output: SequenceTypingOutput) -> None:
        """
        Exports the output of the workflow.
        :param output: Output
        :return: None
        """
        self._report.add_html_object(output.report_section)
        output.report_section.copy_files(self._report.output_dir)

        # Add log files
        dir_logs = os.path.join(self._report.output_dir, 'logs')
        if not os.path.isdir(dir_logs):
            os.makedirs(dir_logs)
        shutil.copyfile(output.log_file, os.path.join(dir_logs, 'log_gene_detection.txt'))
        for key, path in self._helper.logs.items():
            shutil.copyfile(path, os.path.join(dir_logs, f'log_{key}.txt'))

        self._report.save()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    sequence_typing = MainSequenceTyping()
    sequence_typing.run()
