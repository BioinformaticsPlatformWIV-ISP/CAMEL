#!/usr/bin/env python
import argparse
import logging
from pathlib import Path
from typing import Optional

import os

from camel.app.camel import Camel
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.mainscripthelper import MainScriptHelper
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.blast.blastn import Blastn
from camel.app.tools.spatyping.spatyping import SpaTyping
from camel.app.tools.spatyping.spatypingreporter import SpaTypingReporter
from camel.resources import CSS_STYLE


class MainSpaTyping(object):
    """
    This tool is used to run the Spa typing tool.
    """

    def __init__(self, args: Optional[argparse.Namespace] = None):
        """
        Initializes the main script.
        :param args: Arguments, if not set they are removed from the command line
        """
        self._camel = Camel.get_instance()
        self._args = MainSpaTyping._parse_arguments() if args is None else args
        self._sample_name = MainScriptHelper.determine_sample_name(self._args)
        self._helper = MainScriptHelper(self._args.working_dir, self._sample_name)
        self._report = None
        self._db_path = Path(self._args.db_path)

    @staticmethod
    def _parse_arguments() -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        MainScriptHelper.add_common_arguments(argument_parser)
        MainScriptHelper.add_input_files_arguments(argument_parser)
        MainScriptHelper.add_assembly_arguments(argument_parser)
        argument_parser.add_argument('--db-path', help="Path to the database", default='/db/pipelines/saureus')
        return argument_parser.parse_args()

    def run(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._report = self._helper.init_report(
            self._args.output_html, self._args.output_dir, 'Spa typing report', f'<i>spa</i> typing')
        self._helper.export_analysis_info_section(self._report, self._helper.determine_input_files(self._args))
        input_files = self._helper.symlink_input_files(self._args.fasta, self._args.fastq_pe)
        fasta_file = self._helper.get_blast_input(input_files, self._args, self._report)
        blastn_tsv_output = self.__run_blastn(fasta_file)
        spa_typing = self.__run_spa_tying(blastn_tsv_output, fasta_file)
        self.__add_report_output(spa_typing)

    def __init_report(self) -> None:
        """
        Initializes the HTML report.
        :return: None
        """
        self._report = HtmlReport(self._args.output_html, self._args.output_dir)
        if not os.path.isdir(self._args.output_dir):
            os.makedirs(self._args.output_dir)
        self._report.initialize('Spa typing', CSS_STYLE)
        self._report.add_pipeline_header('<i>spa</i> typing')
        self._report.save()

    def __run_blastn(self, fasta_file: ToolIOFile) -> ToolIOFile:
        """
        Runs the BLASTN alignment.
        :param fasta_file: Input FASTA file
        :return: None
        """
        blastn = Blastn(self._camel)
        blastn.add_input_files({
            'DB_BLAST': [ToolIOFile(str(self._db_path / 'profiles.fasta'))],
            'FASTA': [fasta_file]})
        blastn.update_parameters(output_format=SpaTyping.BLASTN_OUTPUT_FORMAT)
        blastn.run(self._args.working_dir)
        return blastn.tool_outputs['TSV'][0]

    def __run_spa_tying(self, tsv_output: ToolIOFile, fasta_file: ToolIOFile) -> SpaTyping:
        """
        Runs the Spa typing on the tabular blast output.
        :param tsv_output: Tabular blast output
        :param fasta_file: FASTA file
        :return: None
        """
        spa_typing = SpaTyping(self._camel)
        spa_typing.add_input_files({
            'TSV': [tsv_output],
            'FASTA': [fasta_file],
            'CSV_profiles': [ToolIOFile(str(self._db_path / 'spatypes.csv'))]
        })
        spa_typing.run(self._args.working_dir)
        return spa_typing

    def __add_report_output(self, spa_typing: SpaTyping) -> None:
        """
        Adds the spa typing output to the report.
        :param spa_typing: Spa typing tool instance
        :return: None
        """
        reporter = SpaTypingReporter(self._camel)
        reporter.add_input_informs({'spa_typing': spa_typing.informs})
        reporter.add_input_files({'VAL_hits': spa_typing.tool_outputs['VAL_hits']})
        reporter.run(self._args.working_dir)
        self._report.add_html_object(reporter.tool_outputs['VAL_HTML'][0].value)
        self._report.save()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    main = MainSpaTyping()
    main.run()
