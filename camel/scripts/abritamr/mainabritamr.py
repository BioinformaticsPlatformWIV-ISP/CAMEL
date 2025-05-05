#!/usr/bin/env python
import argparse
from typing import Optional, Sequence

import yaml

from camel.app.camel import Camel
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.abritamr import CONFIG_DATA, SNAKEFILE_MAIN


class MainAbriTAMR(ReportPipeline):
    """
    Main class to run the AbriTAMR standalone pipeline.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        :return: None
        """
        super().__init__('AbriTAMR standalone', '0.2', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return 'AbriTAMR'

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)

    def __construct_config_file(self, input_files: dict[str, list[dict[str, str]]]) -> str:
        """
        Constructs the configuration file.
        :param input_files: Dictionary with the input files (key can only be FASTA).
        :return: Configuration file
        """

        config_data = self.get_template_data(input_files)
        # Add existing config data
        with open(CONFIG_DATA) as handle_in:
            config_data.update(yaml.safe_load(handle_in.read()))
        config_data['abritamr']['species'] = self._args.species
        config_data['analyses'] = ['abritamr']
        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments (optional)
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(argument_parser)
        argument_parser.add_argument('--species', required=True, choices=[
            'Burkholderia_cepacia', 'Acinetobacter_baumannii', 'Streptococcus_pyogenes',
            'Streptococcus_agalactiae', 'Streptococcus_pneumoniae', 'Enterococcus_faecium',
            'Pseudomonas_aeruginosa', 'Staphylococcus_pseudintermedius', 'Clostridioides_difficile',
            'Klebsiella', 'Neisseria', 'Campylobacter', 'Salmonella', 'Escherichia',
            'Staphylococcus_aureus', 'Burkholderia_pseudomallei', 'Enterococcus_faecalis'])
        return argument_parser.parse_args(args)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainAbriTAMR()
    main.run()
