#!/usr/bin/env python
import argparse
from typing import Optional, Sequence
from pathlib import Path

import pkg_resources

from camel.app.camel import Camel
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils


class MainHybridAssemblyPipeline:
    """
    Main class to run the Hybrid assembly pipeline.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        self._args = MainHybridAssemblyPipeline._parse_arguments(args)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--fastq-pe', type=Path, help="Input Fastq PE file", nargs='+')
        argument_parser.add_argument('--fastq-se', type=Path, help="Input Fastq SE files")
        argument_parser.add_argument('--working-dir', type=Path, default=Path.cwd())
        argument_parser.add_argument('--output', type=Path, default="output.tsv")
        argument_parser.add_argument('--ont-qual', type=str, required=True,
                                     choices=['nano-corr', 'nano-hq', 'nano-raw'], default='nano-corr')
        argument_parser.add_argument('--expected-species', type=str, required=True)
        argument_parser.add_argument('--expected-gc-content', type=str, required=True)
        argument_parser.add_argument('--expected-genome-size', type=str, required=True)
        argument_parser.add_argument('--filtlong-keep-percent', type=int)
        argument_parser.add_argument('--freebayes-ploidy', choices=['GRCh37', 'GRCh38', 'X', 'Y', '1'], default='1')
        argument_parser.add_argument('--freebayes-min-alternate-fraction', type=float, default=0.5)
        argument_parser.add_argument('--freebayes-min-alternate-count', type=int, default=10)
        argument_parser.add_argument('--clair3-haploid-precise', action='store_true', default=None)
        argument_parser.add_argument('--clair3-no-phasing', action='store_true', default=None)
        argument_parser.add_argument('--clair3-include-ctgs', action='store_true', default=None)
        argument_parser.add_argument('--clair3-long-indel', action='store_true', default=None)
        argument_parser.add_argument('--threads', type=int, default=8)
        return argument_parser.parse_args(args)

    def run(self) -> None:
        path_config = self.__create_snakemake_config_data()
        path_snakefile = pkg_resources.resource_filename('camel', 'scripts/hybridassemblypipeline/snakefile/main.smk')
        SnakePipelineUtils.run_snakemake(
            path_snakefile, path_config, [], self._args.working_dir, threads=self._args.threads)

    def __create_snakemake_config_data(self) -> str:
        """
        Creates a Snakemake configuration file.
        :return: Config file data
        """
        config = SnakePipelineUtils.generate_config_file({
            'sample_name': 'Sample', 'working_dir': self._args.working_dir, 'name': 'test_sample',
            'pipeline': {
                'name': 'hybrid assembly pipeline',
                'version': '0.1',
                'title': 'Hybrid Assembly Pipeline'
            },
            'input': {
                'illumina': self._args.fastq_pe,
                'ont': self._args.fastq_se,
            },
            'output': self._args.output,
            'filtlong': {
                'keep_percent': self._args.filtlong_keep_percent if self._args.filtlong_keep_percent is not None else 95
            },
            'assembly': {
                'flye': {
                    'genome_size': self._args.expected_genome_size,
                    self._args.ont_qual.replace('-', '_'): True},
                'min_contig_length': 1000
            },
            'polishing': {
                'medaka': {},
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
            }
        }, self._args.working_dir)
        return config


if __name__ == '__main__':
    Camel.get_instance()
    main = MainHybridAssemblyPipeline()
    main.run()
