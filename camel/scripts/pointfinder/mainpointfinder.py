#!/usr/bin/env python
import argparse
import logging

from camel.app.camel import Camel
from camel.app.components.mainscripthelper import MainScriptHelper
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.pointfinder.pointfinder import PointFinder
from camel.app.tools.pointfinder.pointfinderreporter import PointFinderReporter


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
        MainScriptHelper.add_assembly_arguments(argument_parser)
        MainScriptHelper.add_input_files_arguments(argument_parser)
        argument_parser.add_argument('--species', required=True, choices=[
            'campylobacter', 'enterococcus_faecalis', 'enterococcus_faecium', 'escherichia_coli', 'klebsiella',
            'mycobacterium_tuberculosis', 'neisseria_gonorrhoeae', 'plasmodium_falciparum', 'salmonella'])
        return argument_parser.parse_args()

    def run(self) -> None:
        """
        Runs the tool.
        :return: None
        """
        self._report = self._helper.init_report(
            self._args.output_html, self._args.output_dir, 'PointFinder (local) report', f'PointFinder (local)')
        self._helper.export_analysis_info_section(self._report, self._helper.determine_input_files(self._args))
        input_files = self._helper.symlink_input_files(self._args.fasta, self._args.fastq_pe)
        fasta_file = self._helper.get_blast_input(input_files, self._args, self._report)
        pointfinder = self.__run_pointfinder(fasta_file)
        self.__run_reporter(pointfinder)
        all_informs = self._helper.informs + [pointfinder.informs]
        self._report.add_html_object(SnakePipelineUtils.create_commands_section(all_informs, self._args.working_dir))
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
