#!/usr/bin/env python
import argparse
import logging
from pathlib import Path
from typing import Dict, Any, Optional, Sequence

import shutil

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import variant_calling_clair3
from camel.app.tools.clair3.clair3 import Clair3


class MainCalling(object):
    """
    Class to run samtools variant calling using CAMEL.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        """
        self._args = MainCalling._parse_arguments(args)
        self._camel = Camel()

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--bam', type=Path, help="Input BAM file", required=True)
        argument_parser.add_argument('--reference', type=Path, required=True)
        argument_parser.add_argument('--reference-name')
        argument_parser.add_argument('--output', required=True)
        argument_parser.add_argument('--working-dir', type=Path, default=Path.cwd(), required=True)
        argument_parser.add_argument('--platform', type=str, default='ilmn', required=True)
        argument_parser.add_argument('--model-path', type=Path,
                                     default='/usr/local/bin/lmod/clair3/0.1.12/bin/models/ilmn', required=True)
        argument_parser.add_argument(
            '--output-consensus', help="If specified, the consensus sequence is saved in this file.")
        argument_parser.add_argument('--threads', type=int, default=8)
        argument_parser.add_argument('--haploid-precise', action="store_true", default=False)
        argument_parser.add_argument('--no-phasing', action="store_true", default=False)
        argument_parser.add_argument('--include-ctgs', action="store_true", default=False)
        argument_parser.add_argument('--long-indel', action="store_true", default=False)

        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs the variant calling Snakefile to call the variants.
        :return: None
        """

        clair3 = Clair3(Camel.get_instance())
        clair3.add_input_files({'FASTA': [ToolIOFile(self._args.reference)], 'BAM': [ToolIOFile(self._args.bam)]})
        clair3.update_parameters(platform=self._args.platform, model_path=str(self._args.model_path),
                                 threads=self._args.threads, output_path=str(self._args.output))

        if self._args.haploid_precise:
            clair3.update_parameters(haploid_precise='')
        if self._args.no_phasing:
            clair3.update_parameters(no_phasing='')
        if self._args.include_ctgs:
            clair3.update_parameters(include_ctgs='')
        if self._args.long_indel:
            clair3.update_parameters(long_indel='')

        # clair3.run(self._args.working_dir)
        #     exit()
        #
        # Create config file
        config_data = self.__create_snakemake_config_data()
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

        # Copy input BAM file to the right location
        target_dir = self._args.working_dir / 'variant_calling' / 'read_mapping'
        if not target_dir.exists():
            target_dir.mkdir(parents=True)
        SnakemakeUtils.dump_object([ToolIOFile(self._args.bam)], target_dir / 'bam.io')

        # Run Snakemake to generate output file
        output_path = self._args.working_dir / variant_calling_clair3.OUTPUT_VARIANT_CALLING_UNFILTERED_VCF
        SnakePipelineUtils.run_snakemake(
            variant_calling_clair3.SNAKEFILE_VARIANT_CALLING, config_file, [output_path], self._args.working_dir,
            self._args.threads)

        # # Generate consensus sequence
        # if self._args.output_consensus:
        #     self.__generate_consensus_sequence(self._args.output_consensus, config_data)

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
            'sample_name': 'Sample', 'working_dir': str(self._args.working_dir), 'variant_calling': {
                'platform': self._args.platform,
                'reference': {
                    'name': self._args.reference_name if self._args.reference_name else self._args.reference.name,
                    'path': str(self._args.reference)}
            },
            'model_path': str(self._args.model_path)
        }
        for k in ['haploid_precise', 'no_phasing', 'include_ctgs', 'long_indel']:
            # if (k in self._args) and (vars(self._args)[k] is not False):
            config_data['variant_calling'][k] = vars(self._args)[k]
        return config_data

    def __generate_consensus_sequence(self, output_path: Path, config_data: Dict[str, Any]) -> None:
        """
        Generates the consensus sequence by applying the detected variants to the reference sequence.
        :param output_path: Output path to save the consensus sequence
        :param config_data: Snakemake config data
        :return: None
        """
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir, 'consensus.yml')
        output_path_consensus = self._args.working_dir / variant_calling.OUTPUT_VARIANT_CALLING_CONSENSUS
        SnakePipelineUtils.run_snakemake(
            variant_calling.SNAKEFILE_VARIANT_CALLING, config_file, [output_path_consensus], self._args.working_dir,
            self._args.threads)
        fasta_consensus = SnakemakeUtils.load_object(output_path_consensus)[0].path
        shutil.copyfile(fasta_consensus, output_path)


if __name__ == '__main__':
    main = MainCalling()
    main.run()
