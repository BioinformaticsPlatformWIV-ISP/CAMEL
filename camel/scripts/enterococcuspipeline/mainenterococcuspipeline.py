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
        'kraken2', 'confindr', 'rmlst', 'resfinder', 'ncbi_amr', 'vfdb_core', 'virulencefinder', 'mlst',
        'cgmlst', 'lrefinder', 'plasmidfinder', 'mobsuite', 'bacmet']

    DATA_BY_SPECIES = {
        'faecalis': {
            'cgmlst_db': '/db/sequence_typing/enterococcus_faecalis/cgmlst',
            'full_name': 'Enterococcus faecalis',
            'gc_content': 37.4,
            'genome_size': 2_973_380,
            'mlst_dbs': ['/db/sequence_typing/enterococcus_faecalis/mlst']
            ,
            'pointfinder_db': 'enterococcus_faecalis',
            'quast_fasta': '/db/refgenomes/Enterococcus_faecalis/KB944666.1.fasta',
            'quast_gff': '/db/refgenomes/Enterococcus_faecalis/KB944666.1.gff3'
        },
        'faecium': {
            'cgmlst_db': '/db/sequence_typing/enterococcus_faecium/cgmlst',
            'full_name': 'Enterococcus faecium',
            'gc_content': 38.1,
            'genome_size': 2_796_178,
            'mlst_dbs': [
                '/db/sequence_typing/enterococcus_faecium/mlst',
                '/db/sequence_typing/enterococcus_faecium/mlst_bezdicek'
            ],
            'pointfinder_db': 'enterococcus_faecium',
            'quast_fasta': '/db/refgenomes/Enterococcus_faecium/CP038996.1.fasta',
            'quast_gff': '/db/refgenomes/Enterococcus_faecium/CP038996.1.gff3'
        }
    }

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Enterococcus pipeline', '0.2', SNAKEFILE_MAIN, args)

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
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)

    def __construct_config_file(self, input_files: List[Dict[str, str]]) -> str:
        """
        Constructs the configuration file.
        :return: Configuration file
        """
        config_data = self.get_template_data('fastq_pe', input_files)
        config_data['analyses'] = [key for key in MainEnterococcusPipeline.CUSTOM_ANALYSES if vars(self._args)[key]]
        with CONFIG_DATA.open() as handle_in:
            mainscriptutils.dict_merge(config_data, yaml.load(handle_in.read().format(
                cgmlst_db=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['cgmlst_db'],
                coverage_max=self._args.cov_max,
                expected_species=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['full_name'],
                export_bam='true' if self._args.report_include_bam else 'false',
                export_fastq='true' if self._args.report_include_fastq else 'false',
                gc_content=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['gc_content'],
                genome_size=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['genome_size'],
                mlst_db=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['mlst_db'],
                pointfinder_db=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['pointfinder_db'],
                qc_typing_scheme='cgmlst' if self._args.cgmlst else 'mlst',
                quast_fasta=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['quast_fasta'],
                quast_gff=MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['quast_gff']
            ), Loader=yaml.SafeLoader))

            # Set the species
            config_data['selected_species'] = MainEnterococcusPipeline.DATA_BY_SPECIES[self._args.species]['full_name']

            # Add the kmer-option
            if self._args.spades_kmers is not None:
                config_data['assembly']['spades']['kmers'] = self._args.spades_kmers
                config_data['plasmidspades'] = {'spades': {'kmers': self._args.spades_kmers}}

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
        parser.add_argument('--spades-kmers',
                            help="Comma separated list of K-mers to use for the SPAdes assembly (if not set they are "
                                 "automatically determined by SPAdes)")
        return parser.parse_args(args)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainEnterococcusPipeline()
    main.run()
