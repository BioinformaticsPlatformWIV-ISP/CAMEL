#!/usr/bin/env python
import argparse
import logging
import shutil
from pathlib import Path
from typing import Dict, Any, Optional, Sequence

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources.snakefile import variant_calling_clair3


class MainCalling(object):
    """
    Class to run clair3 variant calling using CAMEL.
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
        argument_parser.add_argument('--working-dir', type=Path, default=Path.cwd())
        argument_parser.add_argument('--platform', type=str, default='ilmn', required=True,
                                     choices=['ont', 'hifi', 'ilmn'])
        argument_parser.add_argument('--model-path', type=str,
                                     default='/usr/local/bin/lmod/clair3/0.1.12/bin/models/ilmn', required=True)
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
        # Create config file
        config_data = self.__create_snakemake_config_data()
        config_file = SnakePipelineUtils.generate_config_file(config_data, self._args.working_dir)

        # Copy input BAM file to the right location
        target_dir = self._args.working_dir / 'variant_calling' / 'read_mapping'
        if not target_dir.exists():
            target_dir.mkdir(parents=True)
        SnakemakeUtils.dump_object([ToolIOFile(self._args.bam)], target_dir / 'bam.io')

        # Run Snakemake to generate output file
        path_vcf = self._args.working_dir / variant_calling_clair3.OUTPUT_VARIANT_CALLING_UNFILTERED_VCF
        SnakePipelineUtils.run_snakemake(
            variant_calling_clair3.SNAKEFILE_VARIANT_CALLING, config_file, [path_vcf], self._args.working_dir,
            self._args.threads)

        # Copy output
        logging.info("Collecting Snakemake output file")
        output_vcf_path = SnakemakeUtils.load_object(path_vcf)[0].path
        shutil.copyfile(output_vcf_path, self._args.output)

    def __create_snakemake_config_data(self) -> Dict:
        """
        Creates a Snakemake configuration file.
        :return: Config file data
        """
        config_data = {
            'sample_name': Path(self._args.bam).stem, 'working_dir': str(self._args.working_dir), 'variant_calling': {
                'platform': self._args.platform,
                'reference': {
                    'name': self._args.reference_name if self._args.reference_name else self._args.reference.name,
                    'path': str(self._args.reference)},
                'bam': self._args.bam
            },
            'model_path': str(self._args.model_path)
        }
        for k in ['haploid_precise', 'no_phasing', 'include_ctgs', 'long_indel']:
            config_data['variant_calling'][k] = vars(self._args)[k]
        return config_data


if __name__ == '__main__':
    main = MainCalling()
    main.run()
