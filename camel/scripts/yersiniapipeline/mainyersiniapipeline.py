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
    Main class to run the Yersinia pipeline
    """

    CUSTOM_ANALYSES = ['kraken2', 'confindr', 'cgmlst', 'mlst'] #TODO: the others

    DATA_BY_SPECIES = {
        #used this as ref genome: https://www.genome.jp/kegg-bin/show_organism?org=yen
        'enterocolitica': {
            'amrfinder_species': 'Yersinia enterocolitica',
            'cgmlst_db': 'TODO',
            'full_name': 'Yersinia enterocolitica',
            'gc_content': 47,
            'genome_size': 4_683_620, #TODO: change if different refseq
            'mlst_db': 'TODO',
            'pointfinder_db': 'TODO',
            #TODO: add to db
            'quast_fasta': '/db/refgenomes/Yersinia_enterocolitica/TODO.fasta',
            'quast_gff': '/db/refgenomes/Yersinia_enterocolitica/TODO.gff3',
            'resfinder4_species': 'Other'
        },
        #used this as ref genome: https://www.genome.jp/kegg-bin/show_organism?org=yps
        'pseudotuberculosis': {
            'amrfinder_species': 'TODO', #not in list curated organisms
            'cgmlst_db': 'TODO',
            'full_name': 'Yersinia pseudotuberculosis',
            'gc_content': 47.5,
            'genome_size': 4_840_898, #TODO: change if different refseq
            'mlst_db': 'TODO',
            'pointfinder_db': 'TODO',
            #TODO: add to db
            'quast_fasta': '/db/refgenomes/Yersinia_pseudotuberculosis/TODO.fasta',
            'quast_gff': '/db/refgenomes/Yersinia_pseudotuberculosis/TODO.gff3',
            'resfinder4_species': 'Other'
        }
    }

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class
        :param args: Arguments (optional)
        """
        super().__init__('Yersinia Pipeline', '0.1', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output
        :return: Title
        """
        return '<i>Yersinia</i> pipeline'

    def run(self) -> None:
        """
        Runs the pipeline
        :return: None
        """
        input_files = self._symlink_input()
        self._validate_fastq_input()
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)
        self._export_assembly()

    def __construct_config_file(self, input_files: List[Dict[str, str]]) -> str:
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
                    amrfinder_species=MainYersiniaPipeline.DATA_BY_SPECIES[self._args.species]['amrfinder_species'],
                    cgmlst_db=MainYersiniaPipeline.DATA_BY_SPECIES[self._args.species]['cgmlst_db'],
                    coverage_max=self._args.cov_max,
                    expected_species=MainYersiniaPipeline.DATA_BY_SPECIES[self._args.species]['full_name'],
                    export_fastq='true' if self._args.report_include_fastq else 'false',
                    gc_content=MainYersiniaPipeline.DATA_BY_SPECIES[self._args.species]['gc_content'],
                    genome_size=MainYersiniaPipeline.DATA_BY_SPECIES[self._args.species]['genome_size'],
                    mlst_db=MainYersiniaPipeline.DATA_BY_SPECIES[self._args.species]['mlst_db'],
                    pointfinder_db=MainYersiniaPipeline.DATA_BY_SPECIES[self._args.species]['pointfinder_db'],
                    #TODO: add the others
                    qc_typing_scheme='cgmlst' if self._args.cgmlst else 'mlst',
                    quast_fasta=MainYersiniaPipeline.DATA_BY_SPECIES[self._args.species]['quast_fasta'],
                    quast_gff=MainYersiniaPipeline.DATA_BY_SPECIES[self._args.species]['quast_gff'],
                    resfinder4_species=MainYersiniaPipeline.DATA_BY_SPECIES[self._args.species][
                        'resfinder4_species'],
                )))
            #TODO: additional MLST schemes

            #set the species
            config_data['selected_species'] = MainYersiniaPipeline.DATA_BY_SPECIES[self._args.species]['full_name']

            #TODO: set the detection method for cgMLST

        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments
        :param args: Command line arguments
        :return: Parsed arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        for analysis_key in MainYersiniaPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_','-')}", action='store_true')
        parser.add_argument('--species', required=True, choices=['enterocolitica', 'pseudotuberculosis'])
        return parser.parse_args()


if __name__ == '__main__':
    Camel.get_instance()
    main = MainYersiniaPipeline()
    main.run()