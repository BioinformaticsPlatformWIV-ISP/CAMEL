#!/usr/bin/env python
import argparse
from pathlib import Path
from typing import Optional, Sequence, List, Dict

import yaml

from camel.app.components import mainscriptutils
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.mockpipeline import SNAKEFILE_MAIN, CONFIG_DATA


class MainMockPipeline(ReportPipeline):
    """
    Base-class for the mock pipeline.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Mock pipeline', '1.0', SNAKEFILE_MAIN, args)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Command line arguments
        :return: Parsed arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        return parser.parse_args(args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return 'Mock pipeline'

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        self._validate_fastq_input()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(str(config_file))
        self._export_assembly()

    def __construct_config_file(self, input_files: List[Dict[str, str]]) -> Path:
        """
        Constructs the configuration file.
        :return: Configuration file
        """
        config_data = self.get_template_data('fastq_pe', input_files)
        config_data['analyses'] = []
        import pprint
        pprint.pprint(config_data)
        exit(0)
        with open(CONFIG_DATA) as handle_in:
            mainscriptutils.dict_merge(
                config_data, yaml.safe_load(handle_in.read().format(
                    coverage_max=self._args.cov_max
                )))
        return Path(SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir))
