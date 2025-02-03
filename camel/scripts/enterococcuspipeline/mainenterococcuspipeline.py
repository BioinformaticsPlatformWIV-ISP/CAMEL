#!/usr/bin/env python
import argparse
from typing import Optional, List, Dict, Sequence, Any

import yaml

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.loggers import logger
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.enterococcuspipeline import SNAKEFILE_MAIN, CONFIG_DATA


class MainEnterococcusPipeline(ReportPipeline):
    """
    Main class to run the Enterococcus pipeline.
    """

    CUSTOM_ANALYSES = [
        'kraken2', 'confindr', 'rmlst', 'lrefinder', 'amrfinder', 'resfinder4', 'vfdb_core', 'virulencefinder', 'mlst',
        'cgmlst', 'plasmidfinder', 'mob_suite', 'bacmet', 'human_read_scrubbing', 'variant_calling']

    DATA_BY_SPECIES = {
        'faecalis': {
            'amrfinder_species': 'Enterococcus_faecalis',
            'cgmlst_db': '/db/sequence_typing/enterococcus_faecalis/cgmlst',
            'full_name': 'Enterococcus faecalis',
            'gc_content': 37.4,
            'genome_size': 2_973_380,
            'mlst_db': '/db/sequence_typing/enterococcus_faecalis/mlst',
            'quast_fasta': '/db/refgenomes/Enterococcus_faecalis/KB944666.1.fasta',
            'quast_gff': '/db/refgenomes/Enterococcus_faecalis/KB944666.1.gff3',
            'reference_name': 'KB944666.1',
            'reference_url': 'https://www.ncbi.nlm.nih.gov/nuccore/KB944666.1',
            'resfinder4_species': 'Enterococcus faecalis'
        },
        'faecium': {
            'amrfinder_species': 'Enterococcus_faecium',
            'cgmlst_db': '/db/sequence_typing/enterococcus_faecium/cgmlst',
            'full_name': 'Enterococcus faecium',
            'gc_content': 38.1,
            'genome_size': 2_796_178,
            'mlst_db': '/db/sequence_typing/enterococcus_faecium/mlst',
            'quast_fasta': '/db/refgenomes/Enterococcus_faecium/CP038996.1.fasta',
            'quast_gff': '/db/refgenomes/Enterococcus_faecium/CP038996.1.gff3',
            'reference_name': 'CP038996.1',
            'reference_url': 'https://www.ncbi.nlm.nih.gov/nuccore/CP038996.1/',
            'resfinder4_species': 'Enterococcus faecium'
        },
        'spp': {
            'amrfinder_species': None,
            'disabled_assays': ['mlst', 'cgmlst', 'variant_calling'],
            'full_name': 'Enterococcus spp.',
            'gc_content': 37.4,
            'genome_size': 2_973_380,
            'resfinder4_species': None
        }
    }

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Enterococcus pipeline', '1.1', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return '<i>Enterococcus</i> pipeline'

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
        config_data['analyses'] = [key for key in MainEnterococcusPipeline.CUSTOM_ANALYSES if vars(self._args)[key]]
        with CONFIG_DATA.open() as handle_in:
            # Note that values are filled in as strings, to get 'None' values in YAML 'null' needs to be used
            mainscriptutils.dict_merge(config_data, yaml.load(handle_in.read().format(
                amrfinder_species=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['amrfinder_species'],
                cgmlst_db=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species].get('cgmlst_db'),
                coverage_max=self._args.cov_max,
                k2_name=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['full_name']
                    if self._args.species != 'spp' else 'Enterococcus',
                k2_level = 'S' if self._args.species != 'spp' else 'G',
                gc_content=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['gc_content'],
                genome_size=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['genome_size'],
                mlst_db=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species].get('mlst_db'),
                qc_typing_scheme='cgmlst' if self._args.cgmlst else 'rmlst',
                quast_fasta=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species].get('quast_fasta', 'null'),
                quast_gff=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species].get('quast_gff', 'null'),
                reference_name=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species].get('reference_name', 'null'),
                reference_url=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species].get('reference_url', 'null'),
                resfinder4_species=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['resfinder4_species'],
                export_bam='true' if self._args.report_include_bam else 'false',
            ), Loader=yaml.SafeLoader))

        # Additional MLST scheme for E. faecium
        if (self._args.species == 'faecium') and self._args.mlst:
            config_data['analyses'].append('mlst_bezdicek')

        # Disable species-specific assays for generic Enterococcus
        config_data['is_generic'] = self._args.species == 'spp'
        if self._args.species == 'spp':
            self._update_config_for_generic_spp(config_data)

        # Set the species
        config_data['selected_species'] = MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['full_name']

        # Set the detection method for cgMLST
        config_data['sequence_typing']['cgmlst']['detection_method'] = {
            'blast': 'blast', 'srst2': 'blast', 'kma': 'kma'}.get(self._args.detection_method)

        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        for analysis_key in MainEnterococcusPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_', '-')}", action='store_true')
        parser.add_argument('--species', required=True, choices=['faecium', 'faecalis', 'spp'])
        return parser.parse_args(args)

    def _update_config_for_generic_spp(self, config_data: Dict[str, Any]) -> None:
        """
        Updates the config file with specific adaptation for generic enterococcus.
        :param config_data: Configuration data
        :return: None
        """
        # Disable incompatible assays
        disabled_assays = MainEnterococcusPipeline.DATA_BY_SPECIES['spp']['disabled_assays']
        config_data['analyses'] = [a for a in config_data['analyses'] if a not in disabled_assays]
        logger.warning(f"Generic 'Enterococcus' selected as species, disabling assays: {', '.join(disabled_assays)}")

        # Disable species specific AMR detection
        config_data['amrfinder']['species'] = None
        config_data['resfinder4']['species'] = None
        config_data['resfinder4']['point'] = False

        # Change the typing scheme for the QC check (no cgMLST is available)
        logger.warning(f"cgMLST is not available for generic 'Enterococcus', using rMLST for the QC check.")
        config_data['quality_checks']['typing_scheme'] = 'rmlst'


if __name__ == '__main__':
    Camel.get_instance()
    main = MainEnterococcusPipeline()
    main.run()
