#!/usr/bin/env python
import argparse
from typing import Optional, Dict, List, Sequence

import yaml

from camel.app.camel import Camel
from camel.app.components.pipelines.basepipeline import BasePipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.influenzapipeline import SNAKEFILE_MAIN, CONFIG_DATA


class MainInfluenzaPipeline(BasePipeline):
    """
    Main class to run the Mycobacterium pipeline.
    """

    CUSTOM_ANALYSES = []

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Influenza pipeline', '0.1', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return 'Influenza pipeline'

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        config_file = self.construct_config_file(input_files)
        self._run_snakemake_main(config_file)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments (optional)
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        BasePipeline.add_common_arguments(parser)
        for analysis_key in MainInfluenzaPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_', '-')}", action='store_true')
        parser.add_argument('--subtype', choices=['A', 'B', 'C'], type=str.upper, default='N/A')
        parser.add_argument('--viral-species', choices=['influenza'], type=str.lower, required=True)
        parser.add_argument('--deconseq-dbs')
        parser.add_argument('--deconseq-sequential', action='store_true')
        parser.add_argument('--deconseq-retain-dbs')
        return parser.parse_args(args)

    def construct_config_file(self, input_files: List[Dict[str, str]]) -> str:
        """
        Constructs the configuration file.
        :return: Configuration file
        """
        config_data = self.get_template_data('fastq_pe', input_files)
        print(config_data)
        config_data.pop('detection_method')
        config_data['analyses'] = [key for key in MainInfluenzaPipeline.CUSTOM_ANALYSES if vars(self._args)[key]]
        with open(CONFIG_DATA) as handle_in:
            config_data.update(yaml.load(handle_in.read().format(
                export_fastq='true' if self._args.report_include_fastq else 'false',
                viral_species=self._args.viral_species,
                subtype=self._args.subtype if self._args.subtype else None,
                deconseq_dbs=self._args.deconseq_dbs if self._args.deconseq_dbs else False,
                deconseq_sequential=self._args.deconseq_sequential,
                deconseq_retain=self._args.deconseq_retain_dbs if self._args.deconseq_retain_dbs else False
            ), Loader=yaml.SafeLoader))
        config_data['quality_checks']['expected_gc_content'] = config_data['gc_content'][self._args.viral_species][self._args.subtype]
        return SnakePipelineUtils.generate_config_file(config_data, self._working_dir)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainInfluenzaPipeline()
    main.run()
