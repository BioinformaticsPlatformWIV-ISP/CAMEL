#!/usr/bin/env python
import argparse
from pathlib import Path
from typing import List, Dict, Optional, Sequence

import yaml

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.loggers import logger
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.bacilluspipeline import SNAKEFILE_MAIN, CONFIG_DATA


class MainBacillusPipeline(ReportPipeline):
    """
    Main class to run the Bacillus pipeline.
    """

    CUSTOM_ANALYSES = {
        'common': ['rmlst', 'plasmidfinder', 'mobsuite', 'vfdb_core', 'amrfinder', 'kraken2', 'confindr', 'straingst'],
        'cereus': ['btyper', 'mlst_cereus', 'cgmlst_cereus'],
        'subtilis': ['fastani', 'mlst_subtilis', 'gmo']
    }

    DATA_BY_SPECIES = {
        'cereus': {
            'gc_content': 35,
            'genome_size': 5_800_000,
            'full_name': 'Bacillus cereus',
            'ref_species': 'NZ_CP017060.1'
        },
        'subtilis': {
            'gc_content': 43,
            'genome_size': 4_200_000,
            'full_name': 'Bacillus subtilis',
            'ref_species': 'NC_000964.3'
        }
    }

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Bacillus pipeline', '0.1', SNAKEFILE_MAIN, args)
        self._args = MainBacillusPipeline._parse_arguments(args)

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
        self._validate_input_files()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)
        self._export_assembly()

    def __get_qc_typing_scheme(self) -> str:
        """
        Returns the typing scheme used for QC.
        :return: Typing scheme
        """
        if self._args.species == 'cereus':
            return 'cgmlst_cereus' if self._args.cgmlst_cereus is True else 'mlst_cereus'
        else:
            return 'mlst_subtilis'

    def __construct_config_file(self, input_files: Dict[str, List[Dict[str, str]]]) -> str:
        """
        Constructs the configuration file.
        :input_files: Dictionary with the input files (keys can be FASTQ_PE, FASTQ_SE).
        :return: Configuration file
        """
        config_data = self.get_template_data(input_files)

        # Analyses to perform
        config_data['analyses'] = []
        for group, keys in MainBacillusPipeline.CUSTOM_ANALYSES.items():
            for key in keys:
                if not vars(self._args)[key]:
                    continue
                if group != 'common' and group != self._args.species:
                    logger.warning(f"Analysis '{key}' not supported for species '{self._args.species}'")
                    continue
                config_data['analyses'].append(key)

        # Parse template
        with open(CONFIG_DATA) as handle_in:
            mainscriptutils.dict_merge(config_data, yaml.load(handle_in.read().format(
                species=self._args.species,
                qc_typing_scheme=self.__get_qc_typing_scheme(),
                coverage_max=self._args.cov_max,
                export_fastq='true' if self._args.report_include_fastq else 'false',
                export_bam='true' if self._args.report_include_bam else 'false',
                expected_gc_content=MainBacillusPipeline.DATA_BY_SPECIES[self._args.species]['gc_content'],
                genome_size=MainBacillusPipeline.DATA_BY_SPECIES[self._args.species]['genome_size'],
                wildcards_assembly='long_read_assembly',
                ref_fasta=Path(Camel.get_instance().config['db_root'], 'refgenomes', 'Bacillus',
                               f'{MainBacillusPipeline.DATA_BY_SPECIES[self._args.species]["ref_species"]}.fasta'),
                ref_gff=Path(Camel.get_instance().config['db_root'], 'refgenomes', 'Bacillus',
                             f'{MainBacillusPipeline.DATA_BY_SPECIES[self._args.species]["ref_species"]}.gff3')
            ), Loader=yaml.SafeLoader))

        # Nanopore settings
        if self._args.input_type in ['nanopore', 'hybrid']:
            config_data['assembly']['flye'] = {
                'genome_size': MainBacillusPipeline.DATA_BY_SPECIES[self._args.species]['genome_size'],
                **config_data['assembly'].get('flye', {})}

        # Disable KMA for hybrid data
        if self._args.input_type == 'hybrid' and self._args.species == 'subtilis':
            config_data['gene_detection']['gmo'].pop('force_detection_method')
            logger.warning('KMA Gene detection is temporary obsolete for hybrid data - reverting to default method')

        # Illumina settings
        if self._args.library is not None:
            config_data['read_trimming']['adapter'] = self._args.library

        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: command line arguments
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        parser.add_argument('--species', type=str, choices=['cereus', 'subtilis'], required=True)
        for _, keys in MainBacillusPipeline.CUSTOM_ANALYSES.items():
            for analysis_key in keys:
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
