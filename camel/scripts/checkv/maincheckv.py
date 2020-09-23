#!/usr/bin/env python
import argparse
from pathlib import Path
from typing import Tuple, Optional, Sequence, Dict, Any

from camel.app.camel import Camel
from camel.app.components.mainscripthelper import MainScriptHelper
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.checkm.checkm import CheckM
from camel.app.tools.checkm.checkmreporter import CheckMReporter
from camel.app.tools.checkv.checkv import CheckV
from camel.app.tools.checkv.checkvreporter import CheckVReporter


class MainCheckV(object):
    """
    This class contains the main script for the CheckV tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        """
        self._args = MainCheckV._parse_arguments(args)
        self._camel = Camel()

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--fasta', nargs=2, help='FASTA input', required=True)
        argument_parser.add_argument('--working-dir', help='Working directory', default=str(Path('.').absolute()))
        argument_parser.add_argument('--output-html', help='Report output')
        argument_parser.add_argument('--output-dir', help='Output directory')
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the tool.
        :return: None
        """
        # Create report
        input_dict, input_files_str = self.__prepare_input()
        report = MainScriptHelper.init_report(self._args.output_html, self._args.output_dir, 'CheckV', 'CheckV')
        MainScriptHelper.export_analysis_info_section(report, input_files_str)

        # Run CheckV
        checkv = CheckV(Camel.get_instance())
        checkv.add_input_files(input_dict)
        checkv.run(self._args.working_dir)

        # Create output report
        checkv_reporter = CheckVReporter(Camel.get_instance())
        checkv_reporter.add_input_files({key: checkv.tool_outputs[key] for key in checkv.tool_outputs.keys()})
        checkv_reporter.run(self._args.working_dir)
        section = checkv_reporter.tool_outputs['HTML'][0].value
        report.add_html_object(section)
        section.copy_files(report.output_dir)

        # Add citation and command
        report.add_html_object(SnakePipelineUtils.create_commands_section([checkv.informs], self._args.working_dir))
        # report.add_html_object(SnakePipelineUtils.create_citations_section(['Parks_2015-checkm']))
        report.save()

    def __prepare_input(self) -> Tuple[Dict[str, Any], str]:
        """
        Prepares the input for the CheckM tool.
        :return: Input dictionary
        """
        dir_input = Path(self._args.working_dir) / 'input'
        dir_input.mkdir(parents=True, exist_ok=True)
        fasta_path, fasta_name = self._args.fasta
        path_new = dir_input / fasta_name
        path_new.symlink_to(fasta_path)
        return {'FASTA': [ToolIOFile(path_new)]}, fasta_name


if __name__ == '__main__':
    Camel.get_instance()
    main = MainCheckV()
    main.run()
