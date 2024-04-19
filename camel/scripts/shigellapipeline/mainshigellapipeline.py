#!/usr/bin/env python
import argparse
from typing import List, Dict, Optional, Sequence

import yaml

from camel.app.camel import Camel
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.shigellapipeline import SNAKEFILE_MAIN, CONFIG_DATA


class MainShigellaPipeline(ReportPipeline):
    """
    Main class to run the Shigella pipeline.
    """

    CUSTOM_ANALYSES = ['kraken2', 'confindr', 'rmlst', 'shigeifinder', 'shigatyper', 'amrfinder', 'resfinder4',
                       'mlst_pasteur', 'mlst_warwick', 'cgmlst', 'virulencefinder', 'mob_suite', 'vfdb_core', 'mykrobe']

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Shigella pipeline', '1.1', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return '<i>Shigella</i> pipeline'

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        self._validate_input_files()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)
        self._export_assembly()

    def __construct_config_file(self, input_files: Dict[str, List[Dict[str, str]]]) -> str:
        """
        Constructs the configuration file to run Snakemake.
        :param input_files: Dictionary with the input files (keys can be FASTQ_PE, FASTQ_SE)
        :return: Configuration file
        """
        config_data = self.get_template_data(input_files)
        config_data['analyses'] = [key for key in MainShigellaPipeline.CUSTOM_ANALYSES if vars(self._args)[key]]
        with open(CONFIG_DATA) as handle_in:
            config_data.update(yaml.load(handle_in.read().format(
                coverage_max=self._args.cov_max,
                export_fastq='true' if self._args.report_include_fastq else 'false'
            ), Loader=yaml.SafeLoader))
        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        for analysis_key in MainShigellaPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_', '-')}", action='store_true')
        return parser.parse_args(args)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainShigellaPipeline()
    main.run()
