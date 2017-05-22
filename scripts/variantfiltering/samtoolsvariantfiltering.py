import argparse
import gzip

import logging
import os

from app.camel import Camel
from app.command.command import Command
from app.io.tooliofile import ToolIOFile
from app.tools.bcftools.bcftoolsindex import BcftoolsIndex
from app.tools.bcftools.bcftoolsview import BcftoolsView
from app.tools.variantfiltering.depthfilter import DepthFilter
from app.tools.variantfiltering.distancefilter import DistanceFilter
from app.tools.variantfiltering.mappingqualityfilter import MappingQualityFilter
from app.tools.variantfiltering.snpqualityfilter import SnpQualityFilter
from app.tools.variantfiltering.zscorefilter import ZScoreFilter


class SamtoolsVariantFiltering(object):
    """
    Class to run the samtools variant filtering using CAMEL. 
    """

    def __init__(self):
        """
        Initializes the main script.
        """
        self._args = SamtoolsVariantFiltering._parse_arguments()
        self._camel = Camel()
        self._destination_path = '.'

    @staticmethod
    def _parse_arguments():
        """
        Parses the command line arguments.
        :return: Arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--vcf', required=True, help="Input VCF file")
        argument_parser.add_argument('--bam')
        argument_parser.add_argument('--output')
        argument_parser.add_argument('--min-total-depth', default=10, type=int)
        argument_parser.add_argument('--min-forward-depth', default=1, type=int)
        argument_parser.add_argument('--min-reverse-depth', default=1, type=int)
        argument_parser.add_argument('--min-snp-quality', default=25, type=float)
        argument_parser.add_argument('--min-mapping-quality', default=30, type=int)
        argument_parser.add_argument('--min-distance', default=10, type=int)
        argument_parser.add_argument('--keep-best', action='store_true')
        argument_parser.add_argument('--min-zscore', default=1.96, type=float)
        argument_parser.add_argument('--y-mult', default=10, type=float)
        return argument_parser.parse_args()

    def run(self):
        """
        Filter the input VCF file.
        :return: None
        """
        vcf_file = self.__compress_vcf_input()
        if self._args.min_total_depth != 0 or self._args.min_forward_depth != 0 or self._args.min_reverse_depth != 0:
            vcf_file = self.__run_depth_filter(vcf_file)
        if self._args.min_snp_quality != 0:
            vcf_file = self.__run_snp_quality_filter(vcf_file)
        if self._args.min_mapping_quality != 0:
            vcf_file = self.__run_mapping_quality_filter(vcf_file)
        if self._args.min_distance != 0:
            vcf_file = self.__run_distance_filter(vcf_file)
        if self._args.bam:
            vcf_file = self.__run_zscore_filter(vcf_file)
        self.__save_to_output_file(vcf_file)

    def __compress_vcf_input(self):
        """
        Compresses the VCF input.
        :return:
        """
        bcftools_view = BcftoolsView(self._camel)
        bcftools_view.add_input_files({
            'VCF': [ToolIOFile(self._args.vcf)]
        })
        bcftools_view.update_parameters(output_format='VCF', compress_output=True)
        bcftools_view.run(self._destination_path)
        return bcftools_view.tool_outputs['VCF_GZ'][0]

    def __remove_vdb_annotation(self, vcf_file):
        """
        Removes the variant distance bias field in the VCF file. This field can cause problems with bcftools view.
        :return: VCF file with the annotation removed.
        """
        c = Command()
        new_vcf_file = os.path.join(os.path.abspath(self._destination_path), 'vcf_file_vdb_removed.vcf.gz')
        c.command = 'ml bcftools; bcftools annotate --output-type z -x INFO/VDB {} > {}'.format(
            vcf_file.path, new_vcf_file)
        c.run_command(self._destination_path)
        return ToolIOFile(new_vcf_file)

    def __run_depth_filter(self, vcf_file):
        """
        Applies the depth filter.
        :param vcf_file: VCF file
        :return: filtered VCF file
        """
        logging.info("Performing depth filtering.")
        depth_filter = DepthFilter(self._camel)
        depth_filter.add_input_files({'VCF_GZ': [vcf_file]})
        depth_filter.update_parameters(min_depth=self._args.min_total_depth,
                                       min_forward_depth=self._args.min_forward_depth,
                                       min_reverse_depth=self._args.min_reverse_depth)
        depth_filter.run(self._destination_path)
        return depth_filter.tool_outputs['VCF_GZ'][0]

    def __run_snp_quality_filter(self, vcf_file):
        """
        Applies the SNP quality filter.
        :param vcf_file: VCF file
        :return: Filtered VCF file
        """
        logging.info("Performing snp quality filtering")
        snp_quality_filter = SnpQualityFilter(self._camel)
        snp_quality_filter.add_input_files({'VCF_GZ': [vcf_file]})
        snp_quality_filter.update_parameters(min_snp_quality=self._args.min_snp_quality)
        snp_quality_filter.run(self._destination_path)
        return snp_quality_filter.tool_outputs['VCF_GZ'][0]

    def __run_mapping_quality_filter(self, vcf_file):
        """
        Applies the mapping quality filter.
        :param vcf_file: VCF file
        :return: Filtered VCF file
        """
        logging.info("Performing mapping quality filtering")
        mapping_quality_filter = MappingQualityFilter(self._camel)
        mapping_quality_filter.add_input_files({'VCF_GZ': [vcf_file]})
        mapping_quality_filter.update_parameters(min_mapping_quality=self._args.min_mapping_quality)
        mapping_quality_filter.run(self._destination_path)
        return mapping_quality_filter.tool_outputs['VCF_GZ'][0]

    def __run_distance_filter(self, vcf_file):
        """
        Applies the distance filter.
        :param vcf_file: VCF file
        :return: Filtered VCF file
        """
        vcf_file = self.__remove_vdb_annotation(vcf_file)
        vcf_file = self.__index_vcf(vcf_file)
        logging.info("Performing distance filtering")
        distance_filter = DistanceFilter(self._camel)
        distance_filter.add_input_files({'VCF_GZ': [vcf_file]})
        distance_filter.update_parameters(min_distance=self._args.min_distance, keep_best=self._args.keep_best)
        distance_filter.run(self._destination_path)
        return distance_filter.tool_outputs['VCF_GZ'][0]

    def __run_zscore_filter(self, vcf_file):
        """
        Applies the Z-Score filter.
        :param vcf_file: VCF file
        :return: Filtered VCF file
        """
        vcf_file = self.__index_vcf(vcf_file)
        zscore_filter = ZScoreFilter(self._camel)
        zscore_filter.add_input_files({'VCF_GZ': [vcf_file], 'BAM': [ToolIOFile(self._args.bam)]})
        zscore_filter.update_parameters(min_zscore=self._args.min_zscore, y_multiplier=self._args.y_mult)
        zscore_filter.run(self._destination_path)
        return zscore_filter.tool_outputs['VCF_GZ'][0]

    def __save_to_output_file(self, vcf_file):
        """
        Saves the content of the VCF file to the output file.
        :param vcf_file: VCF file
        :return: None
        """
        with gzip.open(vcf_file.path) as vcf_content:
            with open(self._args.output, 'w') as handle:
                handle.write(vcf_content.read())

    def __index_vcf(self, vcf_file):
        """
        Indexes a VCF file.
        :param vcf_file: VCF file
        :return: Indexed VCF file
        """
        bcftools_index = BcftoolsIndex(self._camel)
        bcftools_index.add_input_files({'VCF_GZ': [vcf_file]})
        bcftools_index.run(self._destination_path)
        return bcftools_index.tool_outputs['VCF_GZ'][0]

if __name__ == '__main__':
    main = SamtoolsVariantFiltering()
    main.run()
