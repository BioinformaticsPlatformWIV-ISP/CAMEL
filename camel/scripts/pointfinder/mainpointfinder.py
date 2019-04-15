#!/usr/bin/env python
import argparse
import datetime
import logging

import os

from camel.app.camel import Camel
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.mainscripthelper import MainScriptHelper
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.pointfinder.pointfinder import PointFinder
from camel.app.tools.pointfinder.pointfinderreporter import PointFinderReporter
from camel.resources import CSS_STYLE
from camel.resources.javascript import JQUERY_SRC


class MainPointFinder(object):
    """
    This class is used to execute the PointFinder tool.
    """

    def __init__(self, args: argparse.Namespace = None):
        """
        Initializes the main script.
        :param args: Arguments (optional)
        """
        self._args = args if args is not None else MainPointFinder.parse_arguments()
        self._sample_name = MainScriptHelper.determine_sample_name(self._args)
        self._helper = MainScriptHelper(self._args.working_dir, self._sample_name)
        self._report = None

    @staticmethod
    def parse_arguments() -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        MainScriptHelper.add_common_arguments(argument_parser)
        argument_parser.add_argument('--fasta', help="Input FASTA file")
        argument_parser.add_argument('--fasta-name', help="Input FASTA file name")
        argument_parser.add_argument('--fastq-pe', help="Input PE FASTQ files", nargs=2)
        argument_parser.add_argument('--fastq-pe-names', help="Input PE FASTQ file names", nargs=2)
        argument_parser.add_argument('--trim-reads', help="Perform read trimming", action='store_true')
        argument_parser.add_argument('--kmers', help="Kmers to use for assembly")
        argument_parser.add_argument('--report-include-fastq', action='store_true')
        argument_parser.add_argument('--species', required=True, choices=[
            'e.coli', 'gonorrhoeae', 'campylobacter', 'salmonella', 'tuberculosis'])
        return argument_parser.parse_args()

    def run(self) -> None:
        """
        Runs the tool.
        :return: None
        """
        self.__init_report()
        self.__add_analysis_info_section()
        input_files = self._helper.symlink_input_files(self._args.fasta, self._args.fastq_pe)
        fasta_file = self._helper.get_blast_input(input_files, self._args, self._report)
        pointfinder = self.__run_pointfinder(fasta_file)
        self.__run_reporter(pointfinder)
        all_informs = self._helper.informs + [pointfinder.informs]
        self._report.add_html_object(SnakePipelineUtils.create_commands_section(all_informs))
        self._report.save()

    def __init_report(self) -> None:
        """
        Initializes the HTML report
        :return: None
        """
        self._report = HtmlReport(self._args.output_html, self._args.output_dir, include_js=[JQUERY_SRC])
        if not os.path.isdir(self._args.output_dir):
            os.makedirs(self._args.output_dir)
        self._report.initialize('PointFinder (local) report', CSS_STYLE)
        self._report.add_pipeline_header('PointFinder (local)')

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

    def __run_pointfinder(self, fasta_file: ToolIOFile) -> PointFinder:
        """
        Runs the PointFinder tool.
        :param fasta_file: Input FASTA file
        :return: PointFinder tool instance
        """
        camel = Camel()
        pointfinder = PointFinder(camel)
        pointfinder.add_input_files({'FASTA': [fasta_file]})
        pointfinder.update_parameters(database=self._args.species)
        pointfinder.run(self._args.working_dir)
        return pointfinder

    def __run_reporter(self, pointfinder: PointFinder) -> None:
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
        self._report.add_html_object(reporter.tool_outputs['VAL_HTML'][0].value)


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    p = MainPointFinder()
    p.run()
