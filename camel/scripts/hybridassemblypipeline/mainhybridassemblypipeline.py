#!/usr/bin/env python
import argparse
from typing import Optional, Sequence
from pathlib import Path

import pkg_resources
import yaml

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.galaxy.galaxyutils import GalaxyUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.hybridassemblypipeline import CONFIG_DATA


class MainHybridAssemblyPipeline(object):
    """
    Main class to run the hybrid assembly pipeline.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the hybrid assembly pipeline class.
        :param args: arguments to be parsed
        :return: None
        """
        self._args = MainHybridAssemblyPipeline._parse_arguments(args)
        self._sample_name = MainHybridAssemblyPipeline._determine_sample_name(self._args)

    @staticmethod
    def _determine_sample_name(args: argparse.Namespace) -> str:
        if args.sample_name is not None:
            return args.sample_name
        return GalaxyUtils.determine_sample_name_from_fq([fq.name for fq in args.fastq_pe])

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--fastq-pe', type=Path, help='Input Fastq PE file', nargs=2, required=True)
        argument_parser.add_argument('--fastq-pe-names', help='Input Fastq PE file', nargs=2)
        argument_parser.add_argument('--fastq-se', type=Path, help='Input Fastq SE files')
        argument_parser.add_argument('--fastq-se-name', help='Input Fastq SE file names')
        argument_parser.add_argument('--working-dir', type=Path, default=Path.cwd())
        argument_parser.add_argument('--output', type=Path, default='output.tsv')
        argument_parser.add_argument('--output-html', type=Path, default='output.html')
        argument_parser.add_argument('--sample-name', type=str, default='test_sample')
        argument_parser.add_argument('--ont-qual', type=str, required=True,
                                     choices=['nano-corr', 'nano-hq', 'nano-raw'], default='nano-corr')
        argument_parser.add_argument('--expected-genome-size', type=str, required=True)
        argument_parser.add_argument('--filtlong-keep-percent', type=int)
        argument_parser.add_argument('--freebayes-ploidy', choices=['GRCh37', 'GRCh38', 'X', 'Y', '1'], default='1')
        argument_parser.add_argument('--freebayes-min-alternate-fraction', type=float, default=0.5)
        argument_parser.add_argument('--freebayes-min-alternate-count', type=int, default=10)
        argument_parser.add_argument('--clair3-haploid-precise', action='store_true', default=None)
        argument_parser.add_argument('--clair3-no-phasing', action='store_true', default=None)
        argument_parser.add_argument('--clair3-include-ctgs', action='store_true', default=None)
        argument_parser.add_argument('--clair3-long-indel', action='store_true', default=None)
        argument_parser.add_argument('--sniffles-min-svlen', type=int)
        argument_parser.add_argument('--sniffles-mapq', type=int)
        argument_parser.add_argument('--sniffles-min-support', type=int)
        argument_parser.add_argument('--threads', type=int, default=8)
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the hybrid assembly pipeline.
        :return: None
        """
        path_config = self.__create_snakemake_config_data()
        path_snakefile = pkg_resources.resource_filename('camel', 'scripts/hybridassemblypipeline/snakefile/main.smk')
        SnakePipelineUtils.run_snakemake(
            path_snakefile, path_config, [], self._args.working_dir, threads=self._args.threads)

    def __create_snakemake_config_data(self) -> str:
        """
        Creates a Snakemake configuration file.
        :return: Config file data
        """
        config_data = {
            'sample_name': self._sample_name, 'working_dir': self._args.working_dir,
            'pipeline': {
                'name': 'hybrid assembly pipeline',
                'version': '0.1',
                'title': 'Hybrid assembly pipeline'
            },
            'input': {
                'illumina': self._args.fastq_pe,
                'ont': self._args.fastq_se
            },
            'output': self._args.output,
            'output_html': self._args.output_html,
            'filtlong': {
                'keep_percent': self._args.filtlong_keep_percent if self._args.filtlong_keep_percent is not None else 95
            },
            'assembly': {
                'flye': {
                    'genome_size': self._args.expected_genome_size,
                    'nano_corr': True if self._args.ont_qual == 'nano-corr' else False,
                    'nano_hq': True if self._args.ont_qual == 'nano-hq' else False,
                    'nano_raw': True if self._args.ont_qual == 'nano-raw' else False},
                # self._args.ont_qual.replace('-', '_'): True},
                'min_contig_length': 1000
            },
            'polishing': {
                'medaka': {
                    'consensus': {},
                    'stitch': {}
                },
                'polca': {},
                'polypolish': {}
            },
            'read_mapping': {
                'bwa': {},
                'minimap2': {}
            },
            'freebayes': {
                'ploidy': self._args.freebayes_ploidy,
                'min_alternate_fraction': self._args.freebayes_min_alternate_fraction,
                'min_alternate_count': self._args.freebayes_min_alternate_count
            },
            'clair3': {
                'haploid_precise': True if self._args.clair3_haploid_precise is not None else False,
                'no_phasing': True if self._args.clair3_no_phasing is not None else False,
                'include_ctgs': True if self._args.clair3_include_ctgs is not None else False,
                'long_indel': True if self._args.clair3_long_indel is not None else False
            },
            'sniffles': {
                'mapq': self._args.sniffles_mapq if self._args.sniffles_mapq is not None else 25,
                'min_support': self._args.sniffles_min_support if self._args.sniffles_min_support is not None else 'auto',
                'min_svlen': self._args.sniffles_min_svlen if self._args.sniffles_min_svlen is not None else 35
            }
        }
        with open(CONFIG_DATA) as handle_in:
            mainscriptutils.dict_merge(
                config_data, yaml.safe_load(handle_in.read()))
        return SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)


if __name__ == '__main__':
    Camel.get_instance()
    main = MainHybridAssemblyPipeline()
    main.run()
