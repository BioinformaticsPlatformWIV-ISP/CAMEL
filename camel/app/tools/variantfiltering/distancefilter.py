import logging
from pathlib import Path
from typing import List, Tuple

import random
import vcf
# noinspection PyProtectedMember
from vcf.model import _Record as Record
# noinspection PyProtectedMember
from vcf.parser import _Filter as Filter

from camel.app.camel import Camel
from camel.app.command.command import Command
from camel.app.tools.variantfiltering.basefilter import BaseFilter


class DistanceFilter(BaseFilter):
    """
    Filters variants based on distance.

    Note:
        This code does not check for duplicate SNP positions. If this is the case the tool will only keep the last one
        in the VCF file.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Variant Filter: Distance', '0.1', camel)

    @property
    def full_name(self) -> str:
        """
        Returns the full name for this filter.
        :return: Full name
        """
        return 'Distance'

    @property
    def description(self) -> str:
        """
        Returns the description for this filter.
        :return: Description
        """
        if 'keep_best' in self._parameters:
            extra = 'Variant with highest SNP quality is kept'
        else:
            extra = 'Both variants are removed'
        return 'Distance ≥<b>{}</b> between variants ({})'.format(
            self._parameters['min_distance'].value, extra)

    def _apply_filter(self) -> None:
        """
        Applies the filtering on the variants.
        :return: None
        """
        self.__set_seed()
        with open(self._tool_inputs['VCF_GZ'][0].path, 'rb') as handle:
            vcf_reader = vcf.Reader(handle)
            all_variants = list(vcf_reader)
        filtered_positions = self.__get_filtered_positions(vcf_reader, all_variants)
        output_uncompressed = self.__create_output_file(vcf_reader, all_variants, filtered_positions)
        self._command = Command(f'bgzip {output_uncompressed}')
        self._execute_command()

    def __set_seed(self) -> None:
        """
        Sets a random seed.
        :return: None
        """
        if 'seed' in self._parameters:
            seed = self._parameters['seed'].value
        else:
            seed = random.random()
        logging.info("Seed: {}".format(seed))
        self._informs['seed'] = seed
        random.seed(seed)

    def __create_output_file(self, vcf_reader: vcf.Reader, all_variants: List[Record],
                             filtered_positions: List[Tuple[str, int]]) -> Path:
        """
        Creates the output file.
        :return: Output file (uncompressed VCF)
        """
        vcf_reader.filters['distance'] = Filter('distance', 'minimal distance between SNPs (custom)')
        output_uncompressed = self.output_path.parent / self.output_path.name.replace('.gz', '')
        with output_uncompressed.open('w') as handle:
            writer = vcf.Writer(handle, vcf_reader)
            for variant in all_variants:
                if (variant.CHROM, variant.POS) in filtered_positions:
                    if 'soft_filter' in self._parameters:
                        variant.FILTER = 'distance'
                    else:
                        continue
                writer.write_record(variant)
            writer.close()
        return output_uncompressed

    def __get_filtered_positions(self, vcf_reader: vcf.Reader, variants: List[Record]) -> List[Tuple[str, int]]:
        """
        Returns a list with positions that should be filtered.
        :param vcf_reader: VCF reader
        :param variants: List of variants
        :return: List of positions that are removed by the filter
        """
        removed_positions = []
        interval_size = int(self._parameters['min_distance'].value)
        for contig_name, contig in vcf_reader.contigs.items():
            variants_by_pos = {v.POS: v for v in variants if v.CHROM == contig_name and (
                    v.FILTER is None) or len(v.FILTER) == 0}

            for i in range(0, contig.length - interval_size):
                interval = range(i, i + interval_size)
                variants_in_interval = [variants_by_pos[pos] for pos in interval if pos in variants_by_pos]

                # Do nothing if there is no SNP or only one
                if 0 <= len(variants_in_interval) < 2:
                    continue

                # Multiple SNPs - check which ones should be removed
                if 'keep_best' in self._parameters:
                    best_variant = max(variants_in_interval, key=lambda x: x.QUAL)
                    removed_variants = [v for v in variants_in_interval if v is not best_variant]
                else:
                    removed_variants = variants_in_interval
                removed_positions.extend([(v.CHROM, v.POS) for v in removed_variants])
        return list(set(removed_positions))
