#!/usr/bin/env python
import argparse
import random
from typing import Optional, Dict, List, Sequence

import yaml

from camel.app.camel import Camel
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.influenzapipeline import SNAKEFILE_MAIN, CONFIG_DATA


class MainInfluenzaPipeline(ReportPipeline):
    """
    Main class to run the Influenza pipeline.
    """

    CUSTOM_ANALYSES = ['deconseq', 'genometyping']

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Influenza pipeline', '0.1', SNAKEFILE_MAIN, args)

    @property
    def title(self) -> str:
        """
        Returns the title of the pipeline as it appears in the HTML output.
        :return: Title
        """
        return 'Influenza pipeline'

    def run(self) -> None:
        """
        Runs the pipeline.
        :return: None
        """
        input_files = self._symlink_input()
        config_file = self.construct_config_file(input_files)
        self._run_snakemake_main(config_file)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments (optional)
        :return: Arguments
        """
        parser = argparse.ArgumentParser()
        ReportPipeline.add_common_arguments(parser)
        for analysis_key in MainInfluenzaPipeline.CUSTOM_ANALYSES:
            parser.add_argument(f"--{analysis_key.replace('_', '-')}", action='store_true')
        parser.add_argument('--subtype', choices=['A', 'B', 'C'], type=str.upper, default='N/A')
        parser.add_argument('--viral-species', choices=['influenza'], type=str.lower, required=True)
        parser.add_argument('--deconseq-dbs')
        parser.add_argument('--deconseq-sequential', action='store_true')
        parser.add_argument('--deconseq-retain-dbs')

        # Genome typing options
        parser.add_argument('--blastn-genometyping-idcutoff', default=97, type=int)
        parser.add_argument('--genometyping-db', choices=['ncbi', 'avian', 'ecdc'])

        # Consensus sequence extraction
        parser.add_argument('--analysis_type', default='assembly', choices=['assembly', 'alignment'])  # assembly / alignment
        parser.add_argument('--assembler', default='spades', choices=['spades', 'velvetoptimiser', 'megahit'])  # spades / velvetoptimiser / megahit
        parser.add_argument('--aligner', default='bwa', choices=['bwa', 'bowtie2'])  # bwa / bowtie2
        parser.add_argument('--variantcaller', default='unifiedgenotyper', choices=['unifiedgenotyper', 'haplotypecaller'])  # unifiedgenotyper / haplotypecaller

        parser.add_argument('--random-seed', type=int)
        return parser.parse_args(args)

    def construct_config_file(self, input_files: List[Dict[str, str]]) -> str:
        """
        Constructs the configuration file.
        :return: Configuration file
        """
        config_data = self.get_template_data('fastq_pe', input_files)
        config_data.pop('detection_method')
        config_data['analyses'] = [key for key in MainInfluenzaPipeline.CUSTOM_ANALYSES if vars(self._args)[key]]
        with open(CONFIG_DATA) as handle_in:
            config_data.update(yaml.load(handle_in.read().format(
                export_fastq='true' if self._args.report_include_fastq else 'false',
                viral_species=self._args.viral_species,
                subtype=self._args.subtype,
                deconseq_dbs=self._args.deconseq_dbs if self._args.deconseq_dbs else False,
                deconseq_sequential=self._args.deconseq_sequential,
                deconseq_retain=self._args.deconseq_retain_dbs if self._args.deconseq_retain_dbs else False,
                blastn_genometyping_idcutoff=self._args.blastn_genometyping_idcutoff
            ), Loader=yaml.SafeLoader))

        virus_name = f'{self._args.viral_species}_{self._args.subtype}' if self._args.subtype != 'N/A' else self._args.viral_species
        species_info = config_data['species_info'].pop(virus_name)
        config_data['species_info'] = species_info

        config_data['quality_checks']['expected_gc_content'] = config_data['species_info']['gc_content']

        if self._args.genometyping:
            config_data['multi_segment'] = ',' in config_data['species_info']['genome_segments']
            config_data['genometyping_db'] = config_data['species_info']['genome_typing_db'][self._args.genometyping_db]
            config_data['genometyping_db_source'] = self._args.genometyping_db
        else:
            config_data['rule_parameters'].pop('blastn_genometyping')

        config_data['random_seed'] = random.randint(1, 10000000) if not self._args.random_seed else self._args.random_seed

        config_data['aligner'] = self._args.aligner

        return SnakePipelineUtils.generate_config_file(config_data, self._working_dir)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainInfluenzaPipeline()
    main.run()
