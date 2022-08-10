#!/usr/bin/env python
import argparse
from pathlib import Path
from typing import List, Dict, Optional, Sequence

import yaml

from camel.app.camel import Camel
from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.components import mainscriptutils
from camel.scripts.bacilluspipeline import SNAKEFILE_MAIN, CONFIG_DATA


class MainBacillusPipeline(ReportPipeline):
    """
    Main class to run the Bacillus pipeline.
    """

    CUSTOM_ANALYSES = ['kraken', 'btyper', 'mlst', 'cgmlst']

    # Not yet up to date!
    DATA_BY_SPECIES = {
        'cereus': {
            'gc_content': 35,
            'genome_size': 2_796_178,
            'full_name': 'Bacillus cereus',
            'mlst_db': '/db/sequence_typing/bacillus_cereus/mlst',
            'cgmlst_db': '/db/sequence_typing/bacillus_cereus/cgmlst'
        }
    }

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Bacillus pipeline', '1.0', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return '<i>Bacillus</i> pipeline'

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
        config_data['analyses'] = [key for key in MainBacillusPipeline.CUSTOM_ANALYSES if vars(self._args)[key]]
        with open(CONFIG_DATA) as handle_in:
            mainscriptutils.dict_merge(config_data, yaml.load(handle_in.read().format(
                qc_typing_scheme='cgmlst' if self._args.cgmlst else 'mlst',
                coverage_max=self._args.cov_max,
                export_fastq='true' if self._args.report_include_fastq else 'false',
                export_bam='true' if self._args.report_include_bam else 'false',
                expected_species=MainBacillusPipeline.DATA_BY_SPECIES[self._args.species]['full_name'],
                expected_gc_content=MainBacillusPipeline.DATA_BY_SPECIES[self._args.species]['gc_content'],
                genome_size=MainBacillusPipeline.DATA_BY_SPECIES[self._args.species]['genome_size'],
                mlst_db=MainBacillusPipeline.DATA_BY_SPECIES[self._args.species]['mlst_db'],
                cgmlst_db=MainBacillusPipeline.DATA_BY_SPECIES[self._args.species]['cgmlst_db']
            ), Loader=yaml.SafeLoader))

        # Set the species
        config_data['selected_species'] = MainBacillusPipeline.DATA_BY_SPECIES[self._args.species]['full_name']

        # Read trimming
        config_data['read_trimming']['export_fastq'] = 'true' if self._args.report_include_fastq else 'false'
        if self._args.library is not None:
            config_data['read_trimming']['adapter'] = self._args.library

        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        parser.add_argument('--fastq-se', type=Path, help="Input SE FASTQ file")
        parser.add_argument('--fastq-se-name', help="Input SE FASTQ file name")
        parser.add_argument('--species', help="Bacillus species under study")

        for analysis_key in MainBacillusPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_', '-')}", action='store_true')
        return parser.parse_args(args)

    @property
    def sample_name(self) -> str:
        """
        Returns the sample name.
        :return: Sample name
        """
        if self._args.fastq_pe is not None:
            return super().sample_name
        else:
            name = self._args.fastq_se_name if (self._args.fastq_se_name is not None) else self._args.fastq_se
            return FastqUtils.get_sample_name(name, FastqUtils.PATTERN_FQ_SE)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainBacillusPipeline()
    main.run()
