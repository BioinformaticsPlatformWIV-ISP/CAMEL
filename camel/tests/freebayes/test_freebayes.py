import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.components.vcf.vcfutils import VCFUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.freebayes.freebayes import Freebayes


class TestFreebayes(CamelTestSuite):
    """
    Tests freebayes.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('clair3')
    FILE_FASTA = ToolIOFile(test_file_dir / 'bsubtilis.fa')
    FILE_BAM_ILLUMINA = ToolIOFile(test_file_dir / 'bsubtilis_illumina.bam')

    def test_freebayes(self) -> None:
        """
        actually testing Freebayes on illumina sequencing data
        """
        freebayes = Freebayes(self.camel)
        freebayes.add_input_files({'FASTA': [TestFreebayes.FILE_FASTA], 'BAM': [TestFreebayes.FILE_BAM_ILLUMINA]})
        freebayes.run(self.running_dir)
        self.verify_output_files(freebayes, 'VCF')
        self.assertGreater(VCFUtils.count_variants(freebayes.tool_outputs['VCF'][0].path), 0)


if __name__ == '__main__':
    unittest.main()
