#!/usr/bin/env python
import argparse
from pathlib import Path
from typing import List, Dict, Optional, Sequence, Tuple

import yaml

from camel.app.camel import Camel
from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.pipelines.reportpipeline import ReportPipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.app.components import mainscriptutils
from camel.scripts.bacilluspipeline import SNAKEFILE_MAIN, CONFIG_DATA


class MainBacillusPipeline(ReportPipeline):
    """
    Main class to run the Bacillus pipeline.
    """

    CUSTOM_ANALYSES = ['kraken', 'btyper', 'mlst', 'cgmlst', 'amrfinder', 'gmo',
                       'vfdb_core', 'plasmidfinder', 'mobsuite']

    DATA_BY_SPECIES = {
        'cereus': {
            'gc_content': 35,
            'genome_size': 2_796_178,
            'full_name': 'Bacillus cereus',
            'mlst_db': '/db/sequence_typing/bacillus_cereus/mlst',
            'cgmlst_db': '/db/sequence_typing/bacillus_cereus/cgmlst'
        },
        'subtilis': {
            'gc_content': 43,
            'genome_size': 4_134_800,
            'full_name': 'Bacillus subtilis',
            'mlst_db': '/db/sequence_typing/bacillus_subtilis/mlst',
            'cgmlst_db': '/db/sequence_typing/bacillus_subtilis/cgmlst'
        }
    }

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        """
        super().__init__('Bacillus pipeline', '1.0', SNAKEFILE_MAIN, args)
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
        config_file = self.__construct_config_file(input_files)
        self._run_snakemake_main(config_file)

    def _get_fastq_input_links(self) -> List[List[Tuple[Path, str]]]:
        """
        Returns the links to the input FASTQ files.
        :return: Links
        """
        links = []
        if self._args.fastq_pe is not None:
            for read_nb, path in enumerate(self._args.fastq_pe, start=1):
                gzipped = FileSystemHelper.is_gzipped(path)
                links.append([path, f"{self.sample_name}_{read_nb}.fastq{'.gz' if gzipped else ''}"])
        else:
            gzipped = FileSystemHelper.is_gzipped(self._args.fastq_se)
            links.append([self._args.fastq_se, f"{self.sample_name}_1.fastq{'.gz' if gzipped else ''}"])
        return links

    def __construct_config_file(self, input_files: List[Dict[str, str]]) -> str:
        """
        Constructs the configuration file.
        :return: Configuration file
        """
        key_fq_in = 'fastq_pe' if (self._args.fastq_pe is not None) else 'fastq_se'
        config_data = self.get_template_data(key_fq_in, input_files)
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
                cgmlst_db=MainBacillusPipeline.DATA_BY_SPECIES[self._args.species]['cgmlst_db'],
            ), Loader=yaml.SafeLoader))

        # Set studies-specific parameters
        config_data['contamination_check']['level_of_depth'] = 'G'
        config_data['read_type'] = self._args.read_type

        # Nanopore settings
        if self._args.read_type == 'nanopore':
            config_data['assembly']['canu'] = {
                'genome_size': MainBacillusPipeline.DATA_BY_SPECIES[self._args.species]['genome_size'],
                **config_data['assembly'].get('canu', {})}
            config_data['quality_checks']['disabled_checks'] = ['coverage', 'fastqc']

        # Read trimming
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
        mainscriptutils.add_input_files_arguments(parser)
        mainscriptutils.add_common_arguments(parser)
        # Logging
        parser.add_argument(
            '--galaxy-job-id', type=str, help='Job id of the run in galaxy (used for logging')
        parser.add_argument(
            '--log', action='store_true', help="If this flag is set, config file and error logs are kept")
        parser.add_argument(
            '--library', help="Adapter library that was used for the sequencing",
            choices=['NexteraPE', 'TruSeq2', 'TruSeq3'], default='NexteraPE')
        parser.add_argument('--output-tsv', help="Output file for the summary", required=True, type=Path)
        parser.add_argument(
            '--report-include-bam', help="Include the BAM file in the report", action='store_true')
        parser.add_argument(
            '--detection-method', help="Type of allele detection: local alignment (blast), read mapping (srst2)",
            choices=['blast', 'kma', 'srst2'], default='blast')
        parser.add_argument(
            '--cov-max', default=100.0, type=float,
            help='Maximum coverage (datasets with higher estimated coverage will be downsampled to the given value)')
        parser.add_argument('--species', help='Bacillus species under study', choices=['cereus', 'subtilis'])

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
