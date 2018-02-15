import logging
import math

import os
import vcf

from app.command.command import Command
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.tools.variantfiltering.filter import Filter


class ZScoreFilter(Filter):
    """
    Filters variants based on Z-score.
    The Z-score is calculated as:
    Z = (X-Y) / sqrt(X+Y)
    Where X is the number of reads having the most common nucleotide at that position and Y the number of reads
    supporting other nucleotides.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(ZScoreFilter, self).__init__('Variant Filter: Z-score', '0.1', camel)

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        if 'BAM' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No BAM input found")
        super(ZScoreFilter, self)._check_input()

    @staticmethod
    def calculate_zscore(x, y):
        """
        Calculates the Z-score.
        :param x: Number of reads with the most common nucleotide.
        :param y: Number of other reads
        :return: Z-score
        """
        return float(x - y) / math.sqrt(x + y)

    @staticmethod
    def get_actg_counts(pileup_line):
        """
        Get the count for each base
        :param pileup_line: Line of the pileup output.
        :return: A, C, T, G counts
        """
        return [pileup_line.upper().count(base) for base in ('A', 'C', 'T', 'G')]

    def __create_pileup(self):
        """
        Creates a pileup file with every position that covers a variant.
        :return: Path to pileup file
        """
        vcf_reader = vcf.Reader(open(self._tool_inputs['VCF_GZ'][0].path, 'rb'))
        positions_file = os.path.join(self._folder, 'positions.txt')
        with open(positions_file, 'w') as handle:
            for variant in vcf_reader:
                handle.write('{}\t{}\n'.format(variant.CHROM, variant.POS))

        path = os.path.join(self._folder, 'positions.pileup')
        pileup_command = Command()
        pileup_command.command = '{}samtools mpileup --count-orphans --positions {} {} > {}'.format(
            self._build_dependencies(), positions_file, self._tool_inputs['BAM'][0].path, path)
        pileup_command.run_command(self._folder)
        return path

    def __filter_positions(self, pileup_file):
        """
        Filters all variant positions based on the Z-score.
        :param pileup_file: Pileup file
        :return: List with the positions that are kept.
        """
        kept_positions = []
        with open(pileup_file, 'r') as handle:
            content = handle.readlines()
            for line in content:
                chrom = line.split('\t')[0]
                pos = line.split('\t')[1]
                bases = line.strip().split('\t')[4]
                actg_counts = ZScoreFilter.get_actg_counts(bases)
                x = max(actg_counts)
                actg_counts.remove(x)
                y = sum(actg_counts)
                zscore = float(x - y) / math.sqrt(x + y)
                logging.debug('Zscore at position {}:{}= {:.2f} {}'.format(chrom, pos, zscore, actg_counts + [x]))

                keep = True
                if zscore < float(self._parameters['min_zscore'].value):
                    logging.info("{}:{} does not pass Z score test".format(chrom, pos))
                    keep = False

                if 'y_multiplier' in self._parameters:
                    if x < float(self._parameters['y_multiplier'].value) * y:
                        logging.info("{}:{} does not pass Y multiplier test".format(chrom, pos))
                        keep = False

                if keep:
                    kept_positions.append('{}\t{}'.format(chrom, pos))
        return kept_positions

    def _apply_filter(self):
        """
        Applies the filtering on the variants.
        :return: None
        """
        pileup_file = self.__create_pileup()
        kept_positions = self.__filter_positions(pileup_file)
        regions_file = os.path.join(self._folder, 'regions_zscore_filter.txt')
        with open(regions_file, 'w') as handle:
            for position in kept_positions:
                handle.write('{}\n'.format(position))
        self.__build_command(regions_file)
        self._execute_command()

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
            self._command.command += ' --regions-file {}'.format(regions_file)
        else:
            self._command.command += ' --exclude 1'
