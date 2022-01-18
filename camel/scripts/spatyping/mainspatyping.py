#!/usr/bin/env python
import argparse
import logging
from pathlib import Path
from typing import Optional, Sequence

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.workflows.readtype import helper_by_read_type
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.blast.blastn import Blastn
from camel.app.tools.spatyping.spatyping import SpaTyping
from camel.app.tools.spatyping.spatypingreporter import SpaTypingReporter


class MainSpaTyping(object):
    """
    This tool is used to run the Spa typing tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        :param args: Arguments, if not set they are removed from the command line
        """
        self._camel = Camel.get_instance()
        self._args = MainSpaTyping._parse_arguments(args)
        self._sample_name = mainscriptutils.determine_sample_name(self._args)
        self._helper = helper_by_read_type[self._args.read_type](Path(self._args.working_dir), self._sample_name)
        self._db_path = Path(self._args.db_path)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        mainscriptutils.add_common_arguments(argument_parser)
        mainscriptutils.add_input_files_arguments(argument_parser)
        mainscriptutils.add_assembly_arguments(argument_parser)
        argument_parser.add_argument('--db-path', help="Path to the database", default='/db/pipelines/saureus')
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Initialize report
        report = mainscriptutils.init_report(
            Path(self._args.output_html), Path(self._args.output_dir), 'spa typing report', f'<i>spa</i> typing')
        report.add_html_object(mainscriptutils.generate_analysis_info_section(self._args))
        report.save()

        # Run the tools
        fasta_file = self._helper.prepare_fasta_input(report, self._args)
        blastn_tsv_output = self.__run_blastn(fasta_file)
        spa_typing = self.__run_spa_tying(Path(blastn_tsv_output.path), fasta_file)
        self.__add_report_output(spa_typing, report)

    def __run_blastn(self, fasta_file: Path) -> ToolIOFile:
        """
        Runs the BLASTN alignment.
        :param fasta_file: Input FASTA file
        :return: None
        """
        blastn = Blastn(self._camel)
        blastn.add_input_files({
            'DB_BLAST': [ToolIOFile(self._db_path / 'profiles.fasta')],
            'FASTA': [ToolIOFile(fasta_file)]})
        blastn.update_parameters(
            output_format=SpaTyping.BLASTN_OUTPUT_FORMAT,
            num_alignments=100_000,
            task='blastn',
            dust='no'
        )
        blastn.run(Path(self._args.working_dir))
        self._helper.informs.append(blastn.informs)
        return blastn.tool_outputs['TSV'][0]

    def __run_spa_tying(self, tsv_output: Path, fasta_file: Path) -> SpaTyping:
        """
        Runs the Spa typing on the tabular blast output.
        :param tsv_output: Tabular blast output
        :param fasta_file: FASTA file
        :return: None
        """
        spa_typing = SpaTyping(self._camel)
        spa_typing.add_input_files({
            'TSV': [ToolIOFile(tsv_output)],
            'FASTA': [ToolIOFile(fasta_file)],
            'CSV_profiles': [ToolIOFile(self._db_path / 'spatypes.csv')]
        })
        spa_typing.run(self._args.working_dir)
        return spa_typing

    def __add_report_output(self, spa_typing: SpaTyping, report: HtmlReport) -> None:
        """
        Adds the spa typing output to the report.
        :param spa_typing: Spa typing tool instance
        :param report: Report to append information to
        :return: None
        """
        reporter = SpaTypingReporter(self._camel)
        reporter.add_input_informs({'spa_typing': spa_typing.informs})
        reporter.add_input_files({'VAL_hits': spa_typing.tool_outputs['VAL_hits']})
        reporter.run(self._args.working_dir)
        self._helper.export_output_and_commands_section(report, reporter.tool_outputs['VAL_HTML'][0].value)
        report.save()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    main = MainSpaTyping()
    main.run()
