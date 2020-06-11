#!/usr/bin/env python
import argparse
import logging
from typing import Optional, List, Dict, Sequence

import yaml

from camel.app.components.pipelines.basepipeline import BasePipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.staphylococcuspipeline import SNAKEFILE_MAIN, CONFIG_DATA


class MainStaphylococcusPipeline(BasePipeline):
    """
    Main class to run the Staphylococcus pipeline.
    """

    CUSTOM_ANALYSES = [
        'kraken', 'resfinder', 'ncbi_amr', 'pointfinder', 'vfdb_core', 'virulencefinder', 'mlst',
        'cgmlst', 'spa_typing', 'sccmec_typing']

    def __init__(self, args: Optional[argparse.Namespace] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Staphylococcus pipeline', '0.1', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return '<i>Staphylococcus</i> pipeline'

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)

    def __construct_config_file(self, input_files: List[Dict[str, str]]) -> str:
        """
        Constructs the configuration file.
        :return: Configuration file
        """
        config_data = self.get_template_data('fastq_pe', input_files)
        config_data['analyses'] = [key for key in MainStaphylococcusPipeline.CUSTOM_ANALYSES if vars(self._args)[key]]
        with CONFIG_DATA.open() as handle_in:
            config_data.update(yaml.load(handle_in.read().format(
                qc_typing_scheme='cgmlst' if self._args.cgmlst else 'mlst',
                export_fastq='true' if self._args.report_include_fastq else 'false',
                export_bam='true' if self._args.report_include_bam else 'false'
            ), Loader=yaml.SafeLoader))
        return SnakePipelineUtils.generate_config_file(config_data, self._working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        BasePipeline.add_common_arguments(parser)
        for analysis_key in MainStaphylococcusPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_', '-')}", action='store_true')
        return parser.parse_args()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    main = MainStaphylococcusPipeline()
    main.run()
