import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.bcftools.bcftoolscsq import BcftoolsCsq
from camel.app.tools.bcftools.bcftoolsfilter import BcftoolsFilter
from camel.app.tools.bcftools.bcftoolsindex import BcftoolsIndex
from camel.app.tools.bcftools.bcftoolsmpileup import BcftoolsMpileup
from camel.app.tools.bcftools.bcftoolsnorm import BcftoolsNorm


class TestBcftools(CamelTestSuite):
    """
    Tests the bcftools tool suite.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('bcftools')
    FILE_VCF_GZ = ToolIOFile(test_file_dir / 'variants.vcf.gz')
    FILE_BED = ToolIOFile(test_file_dir / 'regions.bed')
    FILE_FASTA = ToolIOFile(test_file_dir / 'reference_h37Rv.fasta')
    FILE_GFF = ToolIOFile(test_file_dir / 'annotation_h37Rv.gff')

    FILE_FASTA_TOY = ToolIOFile(test_file_dir / 'toy.fasta')
    FILE_BAM_TOY = ToolIOFile(test_file_dir / 'toy.bam')

    def test_bcftools_csq(self) -> None:
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
        self.verify_output_files(bcftools_csq, 'VCF')

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
        self.verify_output_files(bcftools_filter, 'VCF')

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

    def test_bcftools_norm(self) -> None:
        """
        Tests the bcftools norm command.
        :return: None
        """
        bcftools_norm = BcftoolsNorm(self.camel)
        bcftools_norm.add_input_files({
            'VCF_GZ': [TestBcftools.FILE_VCF_GZ],
            'FASTA': [TestBcftools.FILE_FASTA]
        })
        bcftools_norm.run(self.running_dir)
        self.verify_output_files(bcftools_norm, 'VCF')
        self.assertIn('total', bcftools_norm.informs)

    def test_bcftools_norm_gz_output(self) -> None:
        """
        Tests the bcftools norm command with compressed output.
        :return: None
        """
        bcftools_norm = BcftoolsNorm(self.camel)
        bcftools_norm.add_input_files({
            'VCF_GZ': [TestBcftools.FILE_VCF_GZ],
            'FASTA': [TestBcftools.FILE_FASTA]
        })
        bcftools_norm.update_parameters(output_format='z', rm_dup='snps')
        bcftools_norm.run(self.running_dir)
        self.verify_output_files(bcftools_norm, 'VCF_GZ')
        self.assertIn('total', bcftools_norm.informs)

    def test_bcftools_index(self) -> None:
        """
        Tests the bcftools index function.
        :return: None
        """
        bcftools_index = BcftoolsIndex(self.camel)
        bcftools_index.add_input_files({
            'VCF_GZ': [TestBcftools.FILE_VCF_GZ]
        })
        bcftools_index.run(self.running_dir)
        self.assertTrue(self.running_dir / f'{TestBcftools.FILE_VCF_GZ.path.name}.csi', "CSI index not found")

    def test_bcftools_mpileup(self) -> None:
        """
        Tests the mpileup function of bcftools.
        :return: None
        """
        samtools_mpileup = BcftoolsMpileup(self.camel)
        samtools_mpileup.add_input_files({'BAM': [TestBcftools.FILE_BAM_TOY], 'FASTA': [TestBcftools.FILE_FASTA_TOY]})
        samtools_mpileup.run(self.running_dir)
        self.verify_output_files(samtools_mpileup, 'VCF')


if __name__ == '__main__':
    unittest.main()
