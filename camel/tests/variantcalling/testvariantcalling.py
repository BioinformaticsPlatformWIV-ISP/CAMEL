import argparse
import unittest

import os
import tempfile

from camel.app.camel import Camel
from camel.app.components.vcf.vcfutils import VCFUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.scripts.variantcalling.samtools.maincalling import MainCalling


class TestVariantCalling(unittest.TestCase):
    """
    Tests the variant calling tool.
    """
    camel = Camel()
    running_dir = None

    # Input files
    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'])
    input_bam_file = ToolIOFile(os.path.join(test_file_dir, 'variant_calling', 'toy', 'toy.bam'))
    input_fasta_ref_file = ToolIOFile(os.path.join(test_file_dir, 'variant_calling', 'toy', 'toy.fasta'))

    def setUp(self) -> None:
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(prefix='camel_', dir=TestVariantCalling.camel.config['temp_dir'])

    def test_variant_calling(self) -> None:
        """
        Tests the variant calling main script.
        :return: None
        """
        output_file_vcf = os.path.join(self.running_dir, 'detected_variants.vcf')
        args = argparse.Namespace(
            bam=TestVariantCalling.input_bam_file.path,
            reference=TestVariantCalling.input_fasta_ref_file.path,
            reference_name=os.path.basename(TestVariantCalling.input_fasta_ref_file.path),
            working_dir=self.running_dir,
            output=output_file_vcf,
            output_consensus=None,
            ploidy=1,
            calling_method='consensus',
            skip_variants=None,
            threads=8
        )
        main_calling = MainCalling(args)
        main_calling.run()
        self.assertGreater(os.path.getsize(output_file_vcf), 0)
        self.assertGreater(VCFUtils.count_variants(output_file_vcf), 0)

    def test_variant_calling_consensus(self) -> None:
        """
        Tests the variant calling main script to generate the consensus sequences.
        :return: None
        """
        output_file_vcf = os.path.join(self.running_dir, 'detected_variants.vcf')
        output_file_fasta = os.path.join(self.running_dir, 'consensus.fasta')
        args = argparse.Namespace(
            bam=TestVariantCalling.input_bam_file.path,
            reference=TestVariantCalling.input_fasta_ref_file.path,
            reference_name=os.path.basename(TestVariantCalling.input_fasta_ref_file.path),
            working_dir=self.running_dir,
            output=output_file_vcf,
            output_consensus=output_file_fasta,
            ploidy=1,
            calling_method='consensus',
            skip_variants=None,
            threads=8
        )
        main_calling = MainCalling(args)
        main_calling.run()
        self.assertGreater(os.path.getsize(output_file_fasta), 0)
