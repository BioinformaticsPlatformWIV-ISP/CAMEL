#!/usr/bin/env python
import argparse
from typing import Optional, List, Dict, Sequence

import yaml

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.enterococcuspipeline import SNAKEFILE_MAIN, CONFIG_DATA


class MainEnterococcusPipeline(ReportPipeline):
    """
    Main class to run the Enterococcus pipeline.
    """

    CUSTOM_ANALYSES = [
        'kraken2', 'confindr', 'rmlst', 'lrefinder', 'amrfinder', 'resfinder4', 'vfdb_core', 'virulencefinder', 'mlst',
        'cgmlst', 'plasmidfinder', 'mob_suite', 'bacmet', 'human_read_scrubbing']

    DATA_BY_SPECIES = {
        'faecalis': {
            'amrfinder_species': 'Enterococcus_faecalis',
            'cgmlst_db': '/db/sequence_typing/enterococcus_faecalis/cgmlst',
            'full_name': 'Enterococcus faecalis',
            'gc_content': 37.4,
            'genome_size': 2_973_380,
            'mlst_db': '/db/sequence_typing/enterococcus_faecalis/mlst',
            'pointfinder_db': 'enterococcus_faecalis',
            'quast_fasta': '/db/refgenomes/Enterococcus_faecalis/KB944666.1.fasta',
            'quast_gff': '/db/refgenomes/Enterococcus_faecalis/KB944666.1.gff3',
            'resfinder4_species': 'Enterococcus faecalis'
        },
        'faecium': {
            'amrfinder_species': 'Enterococcus_faecium',
            'cgmlst_db': '/db/sequence_typing/enterococcus_faecium/cgmlst',
            'full_name': 'Enterococcus faecium',
            'gc_content': 38.1,
            'genome_size': 2_796_178,
            'mlst_db': '/db/sequence_typing/enterococcus_faecium/mlst',
            'pointfinder_db': 'enterococcus_faecium',
            'quast_fasta': '/db/refgenomes/Enterococcus_faecium/CP038996.1.fasta',
            'quast_gff': '/db/refgenomes/Enterococcus_faecium/CP038996.1.gff3',
            'resfinder4_species': 'Enterococcus faecium'
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
            mainscriptutils.dict_merge(config_data, yaml.load(handle_in.read().format(
                amrfinder_species=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['amrfinder_species'],
                cgmlst_db=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['cgmlst_db'],
                coverage_max=self._args.cov_max,
                expected_species=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['full_name'],
                export_fastq='true' if self._args.report_include_fastq else 'false',
                gc_content=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['gc_content'],
                genome_size=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['genome_size'],
                mlst_db=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['mlst_db'],
                pointfinder_db=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['pointfinder_db'],
                qc_typing_scheme='cgmlst' if self._args.cgmlst else 'mlst',
                quast_fasta=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['quast_fasta'],
                quast_gff=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['quast_gff'],
                resfinder4_species=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['resfinder4_species'],
            ), Loader=yaml.SafeLoader))

            # Additional MLST scheme for E. faecium
            if (self._args.species == 'faecium') and self._args.mlst:
                config_data['analyses'].append('mlst_bezdicek')

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
        parser.add_argument('--species', required=True, choices=['faecium', 'faecalis'])
        return parser.parse_args(args)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainEnterococcusPipeline()
    main.run()
