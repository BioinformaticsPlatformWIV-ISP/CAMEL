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


class MainCheckM(object):
    """
    This class contains the main script for the CheckM tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        """
        self._args = MainCheckM._parse_arguments(args)
        self._camel = Camel()

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--fasta', nargs=2, action='append', help='FASTA input', required=True)
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
        report = MainScriptHelper.init_report(self._args.output_html, self._args.output_dir, 'CheckM', 'CheckM')
        MainScriptHelper.export_analysis_info_section(report, input_files_str)

        # Run CheckM
        checkm = CheckM(Camel.get_instance())
        checkm.add_input_files(input_dict)
        checkm.run(self._args.working_dir)

        # Create output report
        checkm_reporter = CheckMReporter(Camel.get_instance())
        checkm_reporter.add_input_informs({'checkm': checkm.informs})
        checkm_reporter.add_input_files({'TSV': checkm.tool_outputs['TSV']})
        checkm_reporter.run(self._args.working_dir)
        report.add_html_object(checkm_reporter.tool_outputs['HTML'][0].value)

        # Add citation and command
        report.add_html_object(SnakePipelineUtils.create_commands_section([checkm.informs], self._args.working_dir))
        report.add_html_object(SnakePipelineUtils.create_citations_section(['Parks_2015-checkm']))
        report.save()

    def __prepare_input(self) -> Tuple[Dict[str, Any], str]:
        """
        Prepares the input for the CheckM tool.
        :return: Input dictionary
        """
        dir_input = Path(self._args.working_dir)
        dir_input.mkdir(parents=True, exist_ok=True)
        input_dict = {'FASTA': []}
        for fasta_file, fasta_name in self._args.fasta:
            path_new = (dir_input / fasta_name)
            path_new.symlink_to(fasta_file)
            input_dict['FASTA'].append(ToolIOFile(path_new))
        input_str = ', '.join([fasta_name for _, fasta_name in self._args.fasta])
        return input_dict, input_str


if __name__ == '__main__':
    Camel.get_instance()
    main = MainCheckM()
    main.run()
