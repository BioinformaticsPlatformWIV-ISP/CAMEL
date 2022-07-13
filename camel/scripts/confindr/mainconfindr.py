#!/usr/bin/env python
import argparse
import json
from pathlib import Path
from typing import Optional, Sequence, Dict, Any

import logging

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.confindr.confindr import ConFindr
from camel.app.tools.confindr.confindrreporter import ConFindrReporter


class MainConFindr(object):
    """
    This class contains the main script for the ConFindr tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        """
        self._args = MainConFindr._parse_arguments(args)
        self._camel = Camel()

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        mainscriptutils.add_input_files_arguments(argument_parser, False)
        mainscriptutils.add_common_arguments(argument_parser)
        argument_parser.add_argument(
            '--output-json', type=Path, help='If specified, ConFindr informs are stored in this file')

        # Parameters
        argument_parser.add_argument('--quality-cutoff', type=int, default=20, help='Base quality cutoff')
        argument_parser.add_argument('--base-cutoff', type=int, default=2, help='Number of bases  cutoff')
        argument_parser.add_argument('--base-percentage-cutoff', type=int, default=5, help='Base percentage cutoff')
        argument_parser.add_argument(
            '--min-matching-hashes', type=int, default=150, help='Minimum number of matching KMA hashes')

        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the tool.
        :return: None
        """
        # Initialize report
        report = mainscriptutils.init_report(
            self._args.output_html, self._args.output_dir, 'ConFindr report', 'ConFindr')
        report.add_html_object(mainscriptutils.generate_analysis_info_section(self._args))
        report.save()

        # Run ConFindr
        input_dict = self.__prepare_input()
        confindr = ConFindr(Camel.get_instance())
        confindr.update_parameters(
            quality_cutoff=self._args.quality_cutoff,
            base_cutoff=self._args.base_cutoff,
            base_fraction_cutoff=self._args.base_percentage_cutoff / 100,
            min_matching_hashes=self._args.min_matching_hashes,
            data_type=self._args.read_type.title(),
            threads=self._args.threads
        )
        confindr.add_input_files(input_dict)
        confindr.run(self._args.working_dir)

        # Save informs to file (if specified)
        if self._args.output_json is not None:
            with self._args.output_json.open('w') as handle:
                json.dump(confindr.informs, handle, indent=2)
                logging.info(f'ConFindr informs saved to {self._args.output_json}')

        # Create output report
        confindr_reporter = ConFindrReporter(Camel.get_instance())
        confindr_reporter.add_input_informs({'confindr': confindr.informs})
        confindr_reporter.run(self._args.working_dir)
        report.add_html_object(confindr_reporter.tool_outputs['HTML'][0].value)

        # Add citation and command
        report.add_html_object(SnakePipelineUtils.create_commands_section([confindr.informs], self._args.working_dir))
        report.add_html_object(SnakePipelineUtils.create_citations_section(['Low_2019-confindr', 'Jolley_2012-rmlst']))
        report.save()
        logging.info(f'Report saved to: {self._args.output_html}')

    def __prepare_input(self) -> Dict[str, Any]:
        """
        Prepares the input for the confindr tool.
        :return: Input dictionary
        """
        if self._args.fastq_pe is not None:
            is_gzipped = FileSystemHelper.is_gzipped(self._args.fastq_pe[0])
            input_dict = {f"FASTQ{'_GZ' if is_gzipped else ''}_PE": [ToolIOFile(fq) for fq in self._args.fastq_pe]}
            return input_dict
        else:
            is_gzipped = FileSystemHelper.is_gzipped(self._args.fastq_se)
            input_dict = {f"FASTQ{'_GZ' if is_gzipped else ''}_SE": [ToolIOFile(self._args.fastq_se)]}
            return input_dict


if __name__ == '__main__':
    Camel.get_instance()
    main = MainConFindr()
    main.run()
