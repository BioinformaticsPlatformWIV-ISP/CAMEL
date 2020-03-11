from pathlib import Path

import os
import tempfile
import unittest

from camel.app.camel import Camel
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.bcftools.bcftoolscsq import BcftoolsCsq
from camel.app.tools.bcftools.bcftoolsfilter import BcftoolsFilter


class TestBcftools(unittest.TestCase):
    """
    Tests the bcftools tool suite.
    """

    camel = Camel()
    running_dir = None

    test_file_dir = Path(camel.config['testing']['testfiles_dir']) / 'bcftools'
    FILE_VCF_GZ = ToolIOFile(Path(test_file_dir) / 'variants.vcf.gz')
    FILE_BED = ToolIOFile(Path(test_file_dir) / 'regions.bed')
    FILE_FASTA = ToolIOFile(Path(test_file_dir) / 'reference_h37Rv.fasta')
    FILE_GFF = ToolIOFile(Path(test_file_dir) / 'annotation_h37Rv.gff')

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(None, 'camel_', TestBcftools.camel.config['temp_dir'])

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

    def test_bcftools_filter(self) -> None:
        """
        Tests bcftools filter with valid input.
        :return: None
        """
        bcftools_filter = BcftoolsFilter(self.camel)
        bcftools_filter.add_input_files({
            'VCF_GZ': [TestBcftools.FILE_VCF_GZ],
            'BED': [TestBcftools.FILE_BED]
        })
        bcftools_filter.run(self.running_dir)
        self.assertTrue('VCF' in bcftools_filter.tool_outputs, "No VCF output generated")
        output_file = Path(bcftools_filter.tool_outputs['VCF'][0].path)
        self.assertTrue(output_file.is_file(), "VCF output is not a file")
        self.assertGreater(output_file.stat().st_size, 0, "VCF output file is empty")

    def test_bcftools_filter_invalid_input(self) -> None:
        """
        Tests that bcftools filter fails with invalid input.
        :return: None
        """
        bcftools_filter = BcftoolsFilter(self.camel)
        bcftools_filter.add_input_files({
            'VCF_GZ': [TestBcftools.FILE_FASTA],
            'BED': [TestBcftools.FILE_BED]
        })
        with self.assertRaises(ToolExecutionError):
            bcftools_filter.run(self.running_dir)


if __name__ == '__main__':
    unittest.main()
