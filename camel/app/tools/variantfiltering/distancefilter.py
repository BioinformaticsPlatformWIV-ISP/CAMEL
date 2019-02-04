import logging
from typing import List, Dict, Tuple

import os
import random
import vcf
# noinspection PyProtectedMember
from vcf.model import _Record as Record

from camel.app.camel import Camel
from camel.app.tools.variantfiltering.filter import Filter


class DistanceFilter(Filter):
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
        super(DistanceFilter, self).__init__('Variant Filter: Distance', '0.1', camel)

    @property
    def full_name(self) -> str:
        """
        Returns the full name for this filter.
        :return: Full name
        """
        return 'Distance'

    def _apply_filter(self) -> None:
        """
        Applies the filtering on the variants.
        :return: None
        """
        seed = random.random()
        logging.info("Seed: {}".format(seed))
        self._informs['seed'] = seed
        random.seed(seed)
        regions_filename = self.__get_regions_file()
        self.__build_command(regions_filename)
        self._execute_command()

    @staticmethod
    def __generate_variant_index(variants: List[Record], chrom: str) -> Dict[int, Record]:
        """
        Generates an index for the variants.
        :param variants: Variants
        :param chrom: Chromosome
        :return: Index ({Position: Variant})
        """
        index = {}
        for variant in [v for v in variants if v.CHROM == chrom]:
            index[variant.POS] = variant
        return index

    def __save_regions(self, variants: List[Tuple[str, int]]) -> str:
        """
        Saves the regions in a text file.
        :return: File path
        """
        regions_filename = os.path.join(self._folder, 'removed_regions_distance_filter.txt')
        with open(regions_filename, 'w') as handle:
            for chrom, pos in variants:
                handle.write('{}\t{}'.format(chrom, pos))
                handle.write('\n')
        return regions_filename

    def __get_regions_file(self) -> str:
        """
        Returns the file containing the regions that are kept.
        :return: Regions file path
        """
        with open(self._tool_inputs['VCF_GZ'][0].path, 'rb') as handle:
            vcf_reader = vcf.Reader(handle)
            contigs = vcf_reader.contigs
            variants = list(vcf_reader)

        kept_positions = []
        interval_size = int(self._parameters['min_distance'].value)
        for contig_name, contig in contigs.items():
            variants_by_pos = DistanceFilter.__generate_variant_index(variants, contig_name)
            for i in range(0, contig.length - interval_size):
                interval = range(i, i + interval_size)
                snps_in_interval = [variants_by_pos[pos] for pos in interval if pos in variants_by_pos]

                # Do nothing if there is no SNP or only one
                if 0 <= len(snps_in_interval) < 2:
                    continue

                # Multiple SNPs - check which ones should be removed
                if 'keep_best' in self._parameters:
                    best_snp = max(snps_in_interval, key=lambda x: x.QUAL)
                    removed_snps = [s for s in snps_in_interval if s is not best_snp]
                else:
                    removed_snps = snps_in_interval

                # Remove SNPs
                for removed_snp in removed_snps:
                    # logging.debug("removing SNP {}:{}".format(removed_snp.CHROM, removed_snp.POS))
                    variants_by_pos.pop(removed_snp.POS)

            # Add the kept SNPs to the regions file
            for pos in variants_by_pos.keys():
                kept_positions.append((contig_name, pos))

        return self.__save_regions(kept_positions)

    def __build_command(self, regions_file: str) -> None:
        """
        Builds the command for this tool.
        :param regions_file: File with the included regions
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            self._tool_inputs['VCF_GZ'][0].path,
            '--output-type z',
            '--output {}'.format(self.output_path)
        ])
        if os.path.getsize(regions_file) != 0:
            self._command.command += ' --targets-file {}'.format(regions_file)
        else:
            self._command.command += ' --exclude 1'
