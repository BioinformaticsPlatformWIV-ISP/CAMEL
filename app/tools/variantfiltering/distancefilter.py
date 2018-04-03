import logging

import os
import random
import vcf

from app.tools.variantfiltering.filter import Filter


class DistanceFilter(Filter):
    """
    Filters variants based on distance.
    
    Note:
        This code does not check for duplicate SNP positions. If this is the case the tool will only keep the last one
        in the VCF file.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(DistanceFilter, self).__init__('Variant Filter: Distance', '0.1', camel)

    def _apply_filter(self):
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
    def __generate_variant_index(variants, chrom):
        """
        Generates an index for the variants.
        :param variants: Variants
        :return: Index ({Position: Variant})
        """
        index = {}
        for variant in [v for v in variants if v.CHROM == chrom]:
            index[variant.POS] = variant
        return index

    def __save_regions(self, variants):
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

    def __get_regions_file(self):
        """
        Returns the file containing the regions that are kept.
        :return: Regions file path
        """
        vcf_reader = vcf.Reader(open(self._tool_inputs['VCF_GZ'][0].path, 'rb'))
        contigs = vcf_reader.contigs

        variants = list(vcf_reader)
        kept_positions = []
        for contig_name, contig in contigs.items():
            index = DistanceFilter.__generate_variant_index(variants, contig_name)
            logging.debug("Checking contig: {}".format(contig_name))
            positions = list(range(0, int(self._parameters['min_distance'].value)))
            while positions[-1] < contig.length:
                snps = [index[pos] for pos in positions if pos in index]
                if len(snps) > 1:
                    random.shuffle(snps)
                    if 'keep_best' in self._parameters:
                        best_snp = max(snps, key=lambda x: x.QUAL)
                        removed_snps = [s for s in snps if s is not best_snp]
                    else:
                        removed_snps = snps
                    for removed_snp in removed_snps:
                        logging.debug("removing SNP {}:{}".format(removed_snp.CHROM, removed_snp.POS))
                        index.pop(removed_snp.POS)
                positions.pop(0)
                positions.append(positions[-1] + 1)
            for pos in index.keys():
                kept_positions.append((contig_name, pos))
        return self.__save_regions(kept_positions)

    def __build_command(self, regions_file):
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
