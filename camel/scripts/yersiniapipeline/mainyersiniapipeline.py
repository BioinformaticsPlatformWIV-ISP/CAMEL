#!/usr/bin/env python
import argparse
from typing import Optional, List, Dict, Sequence

import yaml

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.yersiniapipeline import SNAKEFILE_MAIN, CONFIG_DATA


class MainYersiniaPipeline(ReportPipeline):
    """
    Main class to run the Yersinia pipeline.
    """

    CUSTOM_ANALYSES = ['kraken2', 'confindr', 'amrfinder', 'resfinder', 'vfdb_core', 'cgmlst',
                       'mlst', 'mlst_mcnally', 'cgmlst_species', 'cgmlst_yersinia', 'mob_suite']

    DATA_BY_SPECIES = {
        'enterocolitica': {
            'full_name': 'Yersinia enterocolitica',
            'cgmlst_species': '/db/sequence_typing/yersinia/cgmlst_yersinia_enterocolitica'
        },
        'pseudotuberculosis': {
            'full_name': 'Yersinia pseudotuberculosis',
            'cgmlst_species': '/db/sequence_typing/yersinia/cgmlst_yersinia_pseudotuberculosis'
        }
    }

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Yersinia Pipeline', '1.0', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return '<i>Yersinia</i> pipeline'

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
        Constructs the configuration file.
        :param input_files: Dictionary with the input files (keys can be FASTQ_PE, FASTQ_SE).
        :return: Configuration file
        """
        config_data = self.get_template_data(input_files)
        config_data['analyses'] = [key for key in MainYersiniaPipeline.CUSTOM_ANALYSES if vars(self._args)[key]]
        with open(CONFIG_DATA) as handle_in:
            mainscriptutils.dict_merge(
                config_data, yaml.safe_load(handle_in.read().format(
                    qc_typing_scheme='cgmlst' if self._args.cgmlst else 'mlst',
                    export_fastq='true' if self._args.report_include_fastq else 'false',
                    export_bam='true' if self._args.report_include_bam else 'false',
                    coverage_max=self._args.cov_max,
                    cgmlst_species=MainYersiniaPipeline.DATA_BY_SPECIES[self._args.species]['cgmlst_species'],
                )))

            #set the species
            config_data['selected_species'] = MainYersiniaPipeline.DATA_BY_SPECIES[self._args.species]['full_name']

        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Command line arguments
        :return: Parsed arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        for analysis_key in MainYersiniaPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_','-')}", action='store_true')
        parser.add_argument('--species', required=True, choices=['enterocolitica', 'pseudotuberculosis'])
        return parser.parse_args(args)

if __name__ == '__main__':
    Camel.get_instance()
    main = MainYersiniaPipeline()
    main.run()
