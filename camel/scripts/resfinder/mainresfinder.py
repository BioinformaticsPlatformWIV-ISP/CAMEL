#!/usr/bin/env python
import argparse
from pathlib import Path
from typing import Optional, Sequence

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.workflows.readtype import helper_by_read_type
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.resfinder.resfinder import ResFinder
from camel.app.tools.resfinder.resfinderreporter import ResFinderReporter


class MainResFinder(object):
    """
    This class is used to run the main ResFinder local script.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        :param args: Arguments (optional)
        """
        self._args = MainResFinder.parse_arguments(args)
        self._sample_name = mainscriptutils.determine_sample_name(self._args)
        self._helper = helper_by_read_type[self._args.read_type](self._args.working_dir, self._sample_name)

    @staticmethod
    def parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        group = argument_parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--inputfasta', type=Path, help="Input FASTA file")
        group.add_argument('--inputfastq', type=Path, help="Input Fastq files")

        argument_parser.add_argument('--point', action='store_true')
        argument_parser.add_argument('--acquired', action='store_true')

        argument_parser.add_argument('--min_cov', type=float, default=0.6)
        argument_parser.add_argument('--threshold', type=float, default=0.8)

        argument_parser.add_argument('--outputPath', type=Path, help="Output folder")

        argument_parser.add_argument('--species', choices=[
            '"Campylobacter"', '"Campylobacter jejuni"', '"Campylobacter coli"', '"Enterococcus_faecalis"',
            '"Enterococcus faecium"', '"Escherichia coli"', '"Helicobacter pylori"', '"Klebsiella"',
            '"Mycobacterium tuberculosis"', '"Neisseria gonorrhoeae"', '"Plasmodium falciparum"', '"Salmonella"',
            '"Salmonella enterica"', '"Staphylococcus aureus"'])

        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        # Initialize report
        report = mainscriptutils.init_report(
            self._args.output_html, self._args.output_dir, 'ResFinder report', 'ResFinder')
        additional_info = [
            ['Species:', self._args.species if self._args.species is not None else 'Not specified'],
            ['Min % identity:', self._args.threshold if self._args.threshold is not None else 'Curated (default)'],
            ['Min % coverage:', self._args.min_cov],
        ]
        report.add_html_object(mainscriptutils.generate_analysis_info_section(self._args, additional_info))
        report.save()

        # Run tools
        fasta_file = self._helper.prepare_fasta_input(report, self._args)
        resfinder = self.__run_resfinder(fasta_file)
        section = self.__run_reporter(resfinder)
        report.add_html_object(section)
        section.copy_files(report.output_dir)

        # Save report
        all_informs = self._helper.informs + [resfinder.informs]
        report.add_html_object(SnakePipelineUtils.create_commands_section(all_informs, self._args.working_dir))
        report.add_html_object(SnakePipelineUtils.create_citations_section(['Zankari_2012-resfinder']))
        report.save()

    def __run_resfinder(self, fasta_file: Path) -> ResFinder:
        """
        Runs ResFinder
        :param fasta_file: Input FASTA file
        :return: ResFinder tool instance
        """
        camel = Camel()
        resfinder = ResFinder(camel)
        resfinder.add_input_files({'FASTA': [ToolIOFile(fasta_file)]})
        resfinder.update_parameters(output_path=self._args.working_dir, min_cov=0.6, threshold=0.8)
        resfinder.run(self._args.working_dir)
        return resfinder

    def __run_reporter(self, resfinder: ResFinder) -> HtmlReportSection:
        """
        Runs resfinder reporter.
        :param resfinder: ResFinder tool instance
        :return: None
        """
        camel = Camel()
        reporter = ResFinderReporter(camel)
        reporter.add_input_files({'TSV': resfinder.tool_outputs['TSV']})
        reporter.add_input_informs({'resfinder': resfinder.informs})
        reporter.run()
        return reporter.tool_outputs['VAL_HTML'][0].value


if __name__ == '__main__':
    Camel.get_instance()
    resfinder = MainResFinder()
    resfinder.run()
