#!/usr/bin/env python
import argparse
import json
import logging
from typing import Any, Dict, List, Sequence, Optional

import os

from camel.app.components.mainscripthelper import MainScriptHelper
from camel.app.components.workflows.sequencetypingwrapper import SequenceTypingWrapper, SequenceTypingInput, \
    SequenceTypingOutput
from camel.app.io.tooliofile import ToolIOFile


class MainSequenceTyping(object):
    """
    Class to run sequence typing tool, it supports both BLAST+ and SRST2 as detection methods for alleles.
    """

    def __init__(self, args: Optional[Sequence[str]] = None):
        """
        Initializes the main script.
        :param args: (Optional) arguments
        """
        self._args = MainSequenceTyping._parse_arguments(args)
        self._sample_name = MainScriptHelper.determine_sample_name(self._args)
        self._helper = MainScriptHelper(self._args.working_dir, self._sample_name)
        self._report = None

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: (Optional) arguments
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        MainScriptHelper.add_common_arguments(argument_parser)
        MainScriptHelper.add_assembly_arguments(argument_parser)
        MainScriptHelper.add_input_files_arguments(argument_parser)
        argument_parser.add_argument('--scheme-dir', required=True, type=str)
        argument_parser.add_argument('--detection-method', type=str, choices=['blast', 'srst2'], default='blast')
        argument_parser.add_argument('--srst2-max-unaligned-overlap', type=int, default=100)
        return argument_parser.parse_args(args)

    def run(self):
        """
        Runs the workflow.
        :return: None
        """
        self._report = self._helper.init_report(
            self._args.output_html, self._args.output_dir, 'Sequence typing report',
            f'Sequence typing {self._args.detection_method}')
        self._helper.export_analysis_info_section(self._report, self._helper.determine_input_files(self._args))
        input_files = self._helper.symlink_input_files(self._args.fasta, self._args.fastq_pe)
        db_data = MainSequenceTyping.__get_db_metadata(self._args.scheme_dir)
        if self._args.detection_method == 'blast':
            fasta_file = self._helper.get_blast_input(input_files, self._args, self._report)
            output = self.__run_sequence_typing_blast(fasta_file, db_data['name'], self._args.scheme_dir)
        elif self._args.detection_method == 'srst2':
            input_pe = self._helper.get_srst2_input(input_files, self._args, self._report)
            output = self.__run_sequence_typing_srst2(input_pe, db_data['name'], self._args.scheme_dir)
        else:
            raise ValueError(f"Invalid detection method: {self._args.detection_method}")
        self.__export_output(output)

    @staticmethod
    def __get_db_metadata(directory: str) -> Dict[str, Any]:
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
        self._helper.logs['typing'] = output.log_file
        self._helper.export_output_and_commands_section(self._report, output.report_section)


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    sequence_typing = MainSequenceTyping()
    sequence_typing.run()
