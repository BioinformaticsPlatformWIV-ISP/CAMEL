import argparse
import unittest

import os
import tempfile

from camel.app.camel import Camel
from camel.app.components.vcf.vcfutils import VCFUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.variantfiltering.depthfilter import DepthFilter
from camel.app.tools.variantfiltering.distancefilter import DistanceFilter
from camel.app.tools.variantfiltering.mappingqualityfilter import MappingQualityFilter
from camel.app.tools.variantfiltering.snpqualityfilter import SnpQualityFilter
from camel.app.tools.variantfiltering.zscorefilter import ZScoreFilter
from camel.scripts.variantcalling.samtools.mainfiltering import MainFiltering


class TestVariantFiltering(unittest.TestCase):
    """
    Tests the variant filters.
    """

    camel = Camel()
    running_dir = None

    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'], 'variant_calling')

    FILE_VCF_GZ_UNFILTERED = ToolIOFile(os.path.join(test_file_dir, 'unfiltered_variants-myco.vcf.gz'))
    FILE_BAM = ToolIOFile(os.path.join(test_file_dir, 'alignment.bam'))

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(prefix='camel_', dir=TestVariantFiltering.camel.config['temp_dir'])

    def test_depth_filter(self) -> None:
        """
        Tests the depth filter
        :return: None
        """
        depth_filter = DepthFilter(self.camel)
        depth_filter.add_input_files({'VCF_GZ': [TestVariantFiltering.FILE_VCF_GZ_UNFILTERED]})
        depth_filter.update_parameters(min_depth=20, min_forward_depth=2, min_reverse_depth=2)
        depth_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in depth_filter.tool_outputs)
        self.assertGreater(depth_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty")

    def test_mapping_quality_filter(self) -> None:
        """
        Tests the mapping quality filter.
        :return: None
        """
        mq_filter = MappingQualityFilter(self.camel)
        mq_filter.add_input_files({'VCF_GZ': [TestVariantFiltering.FILE_VCF_GZ_UNFILTERED]})
        mq_filter.update_parameters(min_mapping_quality=25)
        mq_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in mq_filter.tool_outputs)
        self.assertGreater(mq_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty")

    def test_snp_quality_filter(self) -> None:
        """
        Tests the SNP quality filter.
        :return: None
        """
        sq_filter = SnpQualityFilter(self.camel)
        sq_filter.add_input_files({'VCF_GZ': [TestVariantFiltering.FILE_VCF_GZ_UNFILTERED]})
        sq_filter.update_parameters(min_snp_quality=50)
        sq_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in sq_filter.tool_outputs)
        self.assertGreater(sq_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty")

    def test_distance_filter(self) -> None:
        """
        Tests the distance filter.
        :return: None
        """
        distance_filter = DistanceFilter(self.camel)
        distance_filter.add_input_files({'VCF_GZ': [TestVariantFiltering.FILE_VCF_GZ_UNFILTERED]})
        distance_filter.update_parameters(min_distance=10, keep_best=True)
        distance_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in distance_filter.tool_outputs)
        self.assertGreater(distance_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty")

    def test_zscore_filter(self) -> None:
        """
        Tests the Z-score filter.
        :return: None
        """
        zscore_filter = ZScoreFilter(self.camel)
        zscore_filter.add_input_files({
            'VCF_GZ': [TestVariantFiltering.FILE_VCF_GZ_UNFILTERED],
            'BAM': [TestVariantFiltering.FILE_BAM]
        })
        zscore_filter.run(self.running_dir)
        self.assertTrue('VCF_GZ' in zscore_filter.tool_outputs)
        self.assertGreater(zscore_filter.tool_outputs['VCF_GZ'][0].size, 0, "Output file is empty")

    def test_variant_filtering_main(self) -> None:
        """
        Tests the main script for the variant filtering.
        :return: None
        """
        number_variants_in = VCFUtils.count_variants(TestVariantFiltering.FILE_VCF_GZ_UNFILTERED.path)
        output_file_vcf = os.path.join(self.running_dir, 'filtered_variants.vcf')
        args = argparse.Namespace(
            vcf=TestVariantFiltering.FILE_VCF_GZ_UNFILTERED.path,
            bam=TestVariantFiltering.FILE_BAM.path,
            working_dir=self.running_dir,
            output_vcf=output_file_vcf,
            min_total_depth=10,
            min_forward_depth=1,
            min_reverse_depth=1,
            min_snp_quality=20,
            min_mapping_quality=25,
            min_distance=10,
            keep_best=True,
            min_zscore=1.96,
            y_mult=4
        )
        main_filtering = MainFiltering(args)
        main_filtering.run()
        self.assertGreater(os.path.getsize(output_file_vcf), 0)
        self.assertLess(VCFUtils.count_variants(output_file_vcf), number_variants_in)


if __name__ == '__main__':
    unittest.main()
