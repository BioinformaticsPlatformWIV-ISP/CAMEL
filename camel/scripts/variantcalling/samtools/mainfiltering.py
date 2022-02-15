#!/usr/bin/env python
import argparse
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import os
import shutil
import yaml

from camel.app.camel import Camel
from camel.app.components.workflows.variantfilteringwrapper import VariantFilteringWrapper


class MainFiltering(object):
    """
    Class to run the samtools variant filtering using CAMEL.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        """
        self._args = MainFiltering._parse_arguments(args)
        self._camel = Camel()

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--vcf', required=True, help="Input VCF file")
        argument_parser.add_argument('--bam')
        argument_parser.add_argument('--output-vcf', required=True)
        argument_parser.add_argument('--output-stats')
        argument_parser.add_argument('--working-dir', default=os.path.abspath('.'))
        argument_parser.add_argument('--min-total-depth', default=10, type=int)
        argument_parser.add_argument('--min-forward-depth', default=1, type=int)
        argument_parser.add_argument('--min-reverse-depth', default=1, type=int)
        argument_parser.add_argument('--min-snp-quality', default=25, type=float)
        argument_parser.add_argument('--min-mapping-quality', default=30, type=int)
        argument_parser.add_argument('--min-distance', default=10, type=int)
        argument_parser.add_argument('--keep-best', action='store_true')
        argument_parser.add_argument('--min-zscore', default=1.96, type=float)
        argument_parser.add_argument('--y-mult', default=10, type=float)
        argument_parser.add_argument('--soft-filter', action='store_true')
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Filter the input VCF file.
        :return: None
        """
        wrapper = VariantFilteringWrapper(Path(self._args.working_dir))
        wrapper.run_workflow(
            Path(self._args.vcf).stem, Path(self._args.vcf), Path(self._args.bam),
            filtering_options=self.__get_filtering_options())
        shutil.copyfile(wrapper.output.vcf_filtered.path, self._args.output_vcf)
        if self._args.output_stats is not None:
            with open(self._args.output_stats, 'w') as handle:
                yaml.dump(wrapper.output.stats, handle)

    def __get_filtering_options(self) -> Dict[str, Any]:
        """
        Returns the dictionary with filtering options.
        :return: Filtering options
        """
        return {
            'soft_filter': self._args.soft_filter,
            'depth': {
                'min_total_depth': self._args.min_total_depth,
                'min_fwd_depth': self._args.min_forward_depth,
                'min_rev_depth': self._args.min_reverse_depth},
            'snp_quality': {
                'min_snp_quality': self._args.min_snp_quality},
            'mapping_quality': {
                'min_mapping_quality': self._args.min_mapping_quality},
            'distance': {
                'min_distance': self._args.min_distance,
                'keep_best': self._args.keep_best},
            'zscore': {
                'min_zscore': self._args.min_zscore,
                'y_multiplier': self._args.y_mult}
        }


if __name__ == '__main__':
    main = MainFiltering()
    main.run()
