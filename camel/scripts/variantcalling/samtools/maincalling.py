import argparse
import logging
from typing import Dict, Any, Optional

import os
import shutil

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import SNAKEFILE_VARIANT_CALLING
from camel.resources.snakefile.variant_calling import OUTPUT_VARIANT_CALLING_UNFILTERED_VCF, \
    OUTPUT_VARIANT_CALLING_CONSENSUS


class MainCalling(object):
    """
    Class to run samtools variant calling using CAMEL.
    """

    def __init__(self, args: Optional[argparse.Namespace] = None) -> None:
        """
        Initializes the main script.
        """
        self._args = MainCalling._parse_arguments() if args is None else args
        self._camel = Camel()

    @staticmethod
    def _parse_arguments() -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        group = argument_parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--bam', help="Input BAM file")
        group.add_argument('--fastq', help="Input Fastq files")
        argument_parser.add_argument('--reference', required=True)
        argument_parser.add_argument('--reference-name', required=True)
        argument_parser.add_argument('--output', required=True)
        argument_parser.add_argument('--output-consensus',
                                     help="If specified, the consensus sequence is saved in this file.")
        argument_parser.add_argument('--working-dir', default=os.path.abspath('.'))
        argument_parser.add_argument('--ploidy', choices=['GRCh37', 'GRCh38', 'X', 'Y', '1'], default='1')
        argument_parser.add_argument('--calling-method', choices=('consensus', 'multiallelic'))
        argument_parser.add_argument('--skip-variants', choices=['snps', 'indels'])
        argument_parser.add_argument('--mutation-rate', type=float)
        argument_parser.add_argument('--minimal-mq', type=int)
        argument_parser.add_argument('--minimal-bq', type=int)
        argument_parser.add_argument('--output-all-sites', action='store_true')
        argument_parser.add_argument('--count-orphans', action='store_true')
        argument_parser.add_argument('--disable-baq', action='store_true')
        argument_parser.add_argument('--threads', type=int, default=8)
        return argument_parser.parse_args()

    def run(self) -> None:
        """
        Runs the variant calling Snakefile to call the variants.
        :return: None
        """
        # Create config file
        config_data = self.__create_snakemake_config_data()

        # Copy input BAM file to the right location
        target_dir = os.path.join(self._args.working_dir, 'variant_calling', 'read_mapping')
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)
        SnakemakeUtils.dump_object([ToolIOFile(self._args.bam)], os.path.join(target_dir, 'bam.io'))

        # Run Snakemake to generate output file
        output_path = os.path.join(self._args.working_dir, OUTPUT_VARIANT_CALLING_UNFILTERED_VCF)
        SnakePipelineUtils.run_snakemake(
            SNAKEFILE_VARIANT_CALLING, config_data, [output_path], self._args.working_dir, self._args.threads)

        # Generate consensus sequence
        if self._args.output_consensus:
            self.__generate_consensus_sequence(self._args.output_consensus, config_data)

        # Copy output
        logging.info("Collecting Snakemake output file")
        output_vcf_path = SnakemakeUtils.load_object(output_path)[0].path
        shutil.copyfile(output_vcf_path, self._args.output)

    def __create_snakemake_config_data(self) -> Dict:
        """
        Creates a Snakemake configuration file.
        :return: Config file data
        """
        config_data = {
            'sample_name': 'Sample', 'working_dir': self._args.working_dir, 'variant_calling': {
                'ploidy': self._args.ploidy,
                'reference': {'fasta': self._args.reference, 'name': self._args.reference,
                              'path': self._args.reference}
            },
            'variant_filtering': {}
        }
        for k in ['calling_method', 'skip_variants', 'mutation_rate', 'minimal_bq', 'minimal_mq', 'count_orphans',
                  'disable_baq']:
            if (k in self._args) and (vars(self._args)[k] is not None):
                config_data['variant_calling'][k] = vars(self._args)[k]
        if ('output_all_sites' in self._args) and self._args.output_all_sites is True:
            config_data['variant_calling']['variants_only'] = False
        return config_data

    def __generate_consensus_sequence(self, output_path: str, config_data: Dict[str, Any]) -> None:
        """
        Generates the consensus sequence by applying the detected variants to the reference sequence.
        :param output_path: Output path to save the consensus sequence
        :param config_data: Snakemake config data
        :return: None
        """
        output_path_consensus = os.path.join(self._args.working_dir, OUTPUT_VARIANT_CALLING_CONSENSUS)
        SnakePipelineUtils.run_snakemake(
            SNAKEFILE_VARIANT_CALLING, config_data, [output_path_consensus], self._args.working_dir, self._args.threads)
        fasta_consensus = SnakemakeUtils.load_object(output_path_consensus)[0].path
        shutil.copyfile(fasta_consensus, output_path)


if __name__ == '__main__':
    main = MainCalling()
    main.run()
