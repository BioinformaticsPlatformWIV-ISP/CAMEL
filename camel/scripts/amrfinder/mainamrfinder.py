#!/usr/bin/env python
import argparse
from pathlib import Path
from typing import Optional, Sequence

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.amrfinder.amrfinder import AMRFinder
from camel.app.tools.amrfinder.amrfinderreporter import AMRFinderReporter


class MainAMRFinder(object):
    """
    Main script for the AMRFinder tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        :param args: Arguments (optional)
        """
        self._args = MainAMRFinder._parse_arguments(args)
        self._sample_name = mainscriptutils.determine_sample_name(self._args)
        self._dir_working = Path(self._args.working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments to parse
        :return: Parsed arguments
        """
        parser = argparse.ArgumentParser()
        mainscriptutils.add_common_arguments(parser)
        parser.add_argument('--fasta', help='Input FASTA file', required=True)
        parser.add_argument('--fasta-name', help="Input FASTA file name (For Galaxy)", type=str)
        parser.add_argument('--db', help='Database', required=True)
        parser.add_argument('--min-cov', help='Minimum target coverage', type=int, default=50)
        parser.add_argument(
            '--min-id', type=int, help='Minimum % identity (defaults to curated cutoffs when not specified)')
        parser.add_argument(
            '--organism', help='Organisms, required to identify point mutations associated with AMR', type=str)
        return parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        # Initialize report
        report = mainscriptutils.init_report(
            Path(self._args.output_html), Path(self._args.output_dir), 'NCBI AMRFinder report', 'NCBI AMRFinder')
        additional_info = [
            ['Organism:', self._args.organism if self._args.organism is not None else 'Not specified'],
            ['Min % identity:', self._args.min_id if self._args.min_id is not None else 'Curated (default)'],
            ['Min % coverage:', self._args.min_cov],
        ]
        report.add_html_object(mainscriptutils.generate_analysis_info_section(self._args, additional_info))
        report.save()

        # Run AMRFinder
        amrfinder = AMRFinder(Camel.get_instance())
        amrfinder.add_input_files({
            'FASTA': [ToolIOFile(Path(self._args.fasta))],
            'DIR': [ToolIODirectory(Path(self._args.db))]
        })
        amrfinder.update_parameters(min_cov=self._args.min_cov / 100.0)
        if self._args.min_id is not None:
            amrfinder.update_parameters(min_ident=self._args.min_id / 100.0)
        if self._args.organism is not None:
            amrfinder.update_parameters(organism=self._args.organism)
        amrfinder.run(self._dir_working)

        # Create output section
        amrfinder_reporter = AMRFinderReporter(Camel.get_instance())
        amrfinder_reporter.add_input_files({'TSV': amrfinder.tool_outputs['TSV']})
        amrfinder_reporter.add_input_informs({'amrfinder': amrfinder.informs})
        amrfinder_reporter.run(self._dir_working)
        report.add_html_object(amrfinder_reporter.tool_outputs['HTML'][0].value)
        amrfinder_reporter.tool_outputs['HTML'][0].value.copy_files(report.output_dir)

        # Add commands and citations
        report.add_html_object(SnakePipelineUtils.create_commands_section([amrfinder.informs], self._dir_working))
        report.add_html_object(SnakePipelineUtils.create_citations_section(['Feldgarden_2019-ndaro']))
        report.save()


if __name__ == '__main__':
    Camel.get_instance()
    main = MainAMRFinder()
    main.run()
