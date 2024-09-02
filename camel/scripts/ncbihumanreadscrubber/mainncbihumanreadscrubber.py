#!/usr/bin/env python
import argparse
from typing import List, Dict, Optional, Sequence

import yaml

from camel.app.camel import Camel
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.ncbihumanreadscrubber import CONFIG_DATA
from camel.scripts.ncbihumanreadscrubber import SNAKEFILE_MAIN


class MainNcbiHumanReadScrubber(ReportPipeline):
    """
    Main class to run the Ncbi human read scrubber tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('NCBI human read scrubber stand alone', '0.2', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return 'NCBI human read scrubber'

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)

    def __construct_config_file(self, input_files: Dict[str, List[Dict[str, str]]]) -> str:
        """
        Constructs the configuration file.
        :return: Configuration file
        """
        config_data = self.get_template_data(input_files)
        config_data['analyses'] = ['human_read_scrubbing']
        config_data['output_removed_reads'] = 'true' if self._args.export_removed_reads else 'false'
        # Add existing config data
        with open(CONFIG_DATA) as handle_in:
            config_data.update(yaml.safe_load(handle_in.read()))
        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(argument_parser)
        argument_parser.add_argument(
            '--export-removed-reads', help="Export the removed reads", action='store_true')
        arguments = argument_parser.parse_args(args)
        # add this input_type line so as to not have to modify the galaxy wrapper for this tool
        arguments.input_type = 'fasta' if arguments.fasta else 'ont' if arguments.fastq_se else arguments.input_type
        return arguments


if __name__ == '__main__':
    Camel.get_instance()
    main = MainNcbiHumanReadScrubber()
    main.run()
