import logging

import math
import os
import vcf

from app.command.command import Command
from app.components.vcf import vcfutils
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.error.invalidparametererror import InvalidParameterError
from app.io.tooliofile import ToolIOFile
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

    def _check_parameters(self):
        """
        Checks the command line parameters.
        :return: None
        """
        if 'min_zscore' not in self._parameters:
            raise InvalidParameterError("Parameter 'min_zscore' not found")
        super(ZScoreFilter, self)._check_parameters()

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

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        nb_of_variants_pre = vcfutils.count_variants(self._tool_inputs['VCF_GZ'][0].path)
        vcf_reader = vcf.Reader(open(self._tool_inputs['VCF_GZ'][0].path, 'r'))
        positions_file = os.path.join(self._folder, 'positions.txt')
        with open(positions_file, 'w') as handle:
            for variant in vcf_reader:
                handle.write('{}\t{}\n'.format(variant.CHROM, variant.POS))

        pileup_file = os.path.join(self._folder, 'positions.pileup')
        pileup_command = Command()
        pileup_command.command = '{}samtools mpileup --count-orphans --positions {} {} > {}'.format(
            self._build_dependencies(), positions_file, self._tool_inputs['BAM'][0].path, pileup_file)
        pileup_command.run_command(self._folder)

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
                    if x < float(self._parameters['y_multiplier'].value)*y:
                        logging.info("{}:{} does not pass Y multiplier test".format(chrom, pos))
                        keep = False

                if keep:
                    kept_positions.append('{}\t{}'.format(chrom, pos))

        regions_file = os.path.join(self._folder, 'regions_zscore_filter.txt')
        with open(regions_file, 'w') as handle:
            for position in kept_positions:
                handle.write('{}\n'.format(position))

        self.__build_command(regions_file)
        self._execute_command()
        output_file = os.path.join(self._folder, self._parameters['output_filename'].value)
        self._tool_outputs['VCF_GZ'] = [ToolIOFile(output_file)]
        nb_of_variants_post = vcfutils.count_variants(output_file)
        logging.info('{}/{} variants passed Z-score filtering'.format(nb_of_variants_post, nb_of_variants_pre))
        self._informs['variants_in'] = nb_of_variants_pre
        self._informs['variants_out'] = nb_of_variants_post

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
            '--output {}'.format(self._parameters['output_filename'].value)
        ])
        if os.path.getsize(regions_file) != 0:
            self._command.command += ' --regions-file {}'.format(regions_file)
        else:
            self._command.command += ' --exclude 1'
