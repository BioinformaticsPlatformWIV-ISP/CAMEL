#!/usr/bin/env python
import argparse
import logging
from pathlib import Path
from typing import Optional, Sequence

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.workflows.readtype import helper_by_read_type
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.pointfinder.pointfinder import PointFinder
from camel.app.tools.pointfinder.pointfinderreporter import PointFinderReporter


class MainPointFinder(object):
    """
    This class is used to execute the PointFinder tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        :param args: Arguments (optional)
        """
        self._args = MainPointFinder.parse_arguments(args)
        self._sample_name = mainscriptutils.determine_sample_name(self._args)
        self._helper = helper_by_read_type[self._args.read_type](Path(self._args.working_dir), self._sample_name)

    @staticmethod
    def parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        mainscriptutils.add_common_arguments(argument_parser)
        mainscriptutils.add_assembly_arguments(argument_parser)
        mainscriptutils.add_input_files_arguments(argument_parser)
        argument_parser.add_argument('--species', required=True, choices=[
            'campylobacter', 'enterococcus_faecalis', 'enterococcus_faecium', 'escherichia_coli', 'klebsiella',
            'mycobacterium_tuberculosis', 'neisseria_gonorrhoeae', 'plasmodium_falciparum', 'salmonella'])
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the tool.
        :return: None
        """
        # Initialize report
        report = mainscriptutils.init_report(
            Path(self._args.output_html), Path(self._args.output_dir), 'PointFinder (local) report',
            'PointFinder (local)')
        report.add_html_object(mainscriptutils.generate_analysis_info_section(self._args))
        report.save()

        # Run the tool and reporter
        fasta_file = self._helper.prepare_fasta_input(report, self._args)
        pointfinder = self.__run_pointfinder(fasta_file)
        section = self.__run_reporter(pointfinder)
        report.add_html_object(section)

        all_informs = self._helper.informs + [pointfinder.informs]
        report.add_html_object(SnakePipelineUtils.create_commands_section(all_informs, self._args.working_dir))
        report.add_html_object(SnakePipelineUtils.create_citations_section(['Zankari_2017-pointfinder']))
        report.save()

    def __run_pointfinder(self, fasta_file: Path) -> PointFinder:
        """
        Runs the PointFinder tool.
        :param fasta_file: Input FASTA file
        :return: PointFinder tool instance
        """
        camel = Camel()
        pointfinder = PointFinder(camel)
        pointfinder.add_input_files({'FASTA': [ToolIOFile(fasta_file)]})
        pointfinder.update_parameters(database=self._args.species)
        pointfinder.run(self._args.working_dir)
        return pointfinder

    def __run_reporter(self, pointfinder: PointFinder) -> HtmlReportSection:
        """
        Runs the PointFinder reporter.
        :param pointfinder: PointFinder tool instance
        :return: None
        """
        camel = Camel()
        reporter = PointFinderReporter(camel)
        reporter.add_input_files({'TSV': pointfinder.tool_outputs['TSV']})
        reporter.add_input_informs({'pointfinder': pointfinder.informs})
        reporter.run()
        return reporter.tool_outputs['VAL_HTML'][0].value


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    p = MainPointFinder()
    p.run()
