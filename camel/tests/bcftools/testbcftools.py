import os
import tempfile
import unittest

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.bcftools.bcftoolscsq import BcftoolsCsq
from camel.app.tools.bcftools.bcftoolsfilter import BcftoolsFilter


class TestBcftools(unittest.TestCase):
    """
    Tests the bcftools.
    """

    camel = Camel()
    running_dir = None

    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'], 'bcftools')

    FILE_VCF_GZ = ToolIOFile(os.path.join(test_file_dir, 'variants.vcf.gz'))
    FILE_BED = ToolIOFile(os.path.join(test_file_dir, 'regions.bed'))
    FILE_FASTA = ToolIOFile(os.path.join(test_file_dir, 'reference_h37Rv.fasta'))
    FILE_GFF = ToolIOFile(os.path.join(test_file_dir, 'annotation_h37Rv.gff'))

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(None, 'camel_', '/scratch/temp')

    def test_bcftools_csq(self):
        """
        Tests bcftools csq.
        :return: None
        """
        bcftools_csq = BcftoolsCsq(self.camel)
        bcftools_csq.add_input_files({
            'VCF_GZ': [TestBcftools.FILE_VCF_GZ],
            'FASTA': [TestBcftools.FILE_FASTA],
            'GFF': [TestBcftools.FILE_GFF]
        })
        bcftools_csq.run(self.running_dir)
        self.assertTrue('VCF' in bcftools_csq.tool_outputs, "No VCF output generated")
        output_file = bcftools_csq.tool_outputs['VCF'][0].path
        self.assertTrue(os.path.isfile(output_file))
        self.assertGreater(os.path.getsize(output_file), 0)

    def test_bcftools_filter(self):
        """
        Tests bcftools filter.
        :return: None
        """
        bcftools_filter = BcftoolsFilter(self.camel)
        bcftools_filter.add_input_files({
            'VCF_GZ': [TestBcftools.FILE_VCF_GZ],
            'BED': [TestBcftools.FILE_BED]
        })
        bcftools_filter.run(self.running_dir)
        self.assertTrue('VCF' in bcftools_filter.tool_outputs, "No VCF output generated")
        output_file = bcftools_filter.tool_outputs['VCF'][0].path
        self.assertTrue(os.path.isfile(output_file))
        self.assertGreater(os.path.getsize(output_file), 0)


if __name__ == '__main__':
    unittest.main()
