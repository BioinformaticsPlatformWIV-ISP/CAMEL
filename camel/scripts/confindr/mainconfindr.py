#!/usr/bin/env python
import argparse
from pathlib import Path
from typing import Tuple, Optional, Sequence, Dict, Any

from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.mainscripthelper import MainScriptHelper
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.tools.confindr.confindr import ConFindr
from camel.app.tools.confindr.confindrreporter import ConFindrReporter


class MainConFindr(object):
    """
    This class contains the main script for the ConFindr tool.
    """

    SNP_MATRIX_FILENAME = 'snp_matrix.fasta'

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

        # PE input
        argument_parser.add_argument('--fastq-pe', nargs=2, help="FASTQ input files")
        argument_parser.add_argument('--fastq-pe-names', nargs=2, help="FASTQ input file names")

        # SE input
        argument_parser.add_argument('--fastq-se', help='FASTQ input file')
        argument_parser.add_argument('--fastq-se-name', help='FASTQ input file name')

        # Directories and output files
        argument_parser.add_argument('--working-dir', help='Working directory', default=str(Path('.').absolute()))
        argument_parser.add_argument('--output-html', help='Report output')
        argument_parser.add_argument('--output-dir', help='Output directory')

        # Parameters
        argument_parser.add_argument('--data-type', type=str, default='illumina', choices=['illumina', 'nanopore'])
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
        # Create report
        input_dict, input_files_str = self.__prepare_input()
        report = MainScriptHelper.init_report(self._args.output_html, self._args.output_dir, 'ConFindr', 'ConFindr')
        MainScriptHelper.export_analysis_info_section(report, input_files_str)

        # Run ConFindr
        confindr = ConFindr(Camel.get_instance())
        confindr.update_parameters(
            quality_cutoff=self._args.quality_cutoff,
            base_cutoff=self._args.base_cutoff,
            base_fraction_cutoff=self._args.base_percentage_cutoff / 100,
            min_matching_hashes=self._args.min_matching_hashes,
            data_type=self._args.data_type.title()
        )
        confindr.add_input_files(input_dict)
        confindr.run(self._args.working_dir)

        # Create output report
        confindr_reporter = ConFindrReporter(Camel.get_instance())
        confindr_reporter.add_input_informs({'confindr': confindr.informs})
        confindr_reporter.run(self._args.working_dir)
        report.add_html_object(confindr_reporter.tool_outputs['HTML'][0].value)

        # Add citation and command
        report.add_html_object(SnakePipelineUtils.create_commands_section([confindr.informs], self._args.working_dir))
        report.add_html_object(SnakePipelineUtils.create_citations_section(['Low_2019-confindr', 'Jolley_2012-rmlst']))
        report.save()

    def __prepare_input(self) -> Tuple[Dict[str, Any], str]:
        """
        Prepares the input for the confindr tool.
        :return: Input dictionary
        """
        if self._args.fastq_pe is not None:
            is_gzipped = FileSystemHelper.is_gzipped(self._args.fastq_pe[0])
            input_dict = {f"FASTQ{'_GZ' if is_gzipped else ''}_PE": [ToolIOFile(fq) for fq in self._args.fastq_pe]}
            names = ', '.join([Path(x).name for x in (
                self._args.fastq_pe_names if self._args.fastq_pe_names else self._args.fastq_pe)])
            return input_dict, names
        else:
            is_gzipped = FileSystemHelper.is_gzipped(self._args.fastq_se)
            input_dict = {f"FASTQ{'_GZ' if is_gzipped else ''}_SE": [ToolIOFile(self._args.fastq_se)]}
            names = self._args.fastq_se_name if self._args.fastq_se_name else Path(self._args.fastq_se).name
            return input_dict, names


if __name__ == '__main__':
    Camel.get_instance()
    main = MainConFindr()
    main.run()
