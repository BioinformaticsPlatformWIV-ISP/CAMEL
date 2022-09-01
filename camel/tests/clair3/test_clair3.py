import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.clair3.clair3 import Clair3


class TestClair3(CamelTestSuite):
    """
    Initializes this testing tool
    """

    test_file_dir = Path('/testdata/camel/btyper/')
    FILE_FASTA_ILLUMINA = ToolIOFile(test_file_dir / 'bacillus_contigs.fasta')
    FILE_BAM_ILLUMINA = ToolIOFile(test_file_dir / 'something.bam')
    FILE_FASTA_ONT = ToolIOFile(test_file_dir / 'bacillus_contigs.fasta')
    FILE_BAM_ONT = ToolIOFile(test_file_dir / 'something.bam')

    def test_clair3_illumina(self) -> None:
        """
        actually testing BTyper
        """
        clair3 = Clair3(self.camel)
        clair3.add_input_files({'FASTA': [TestClair3.FILE_FASTA_ILLUMINA], 'BAM': [TestClair3.FILE_BAM_ILLUMINA]})
        clair3.update_parameters(output_path=self.running_dir, haploid_precise=True, no_phasing=True, include_ctgs=True)
        clair3.run(self.running_dir)
        self.verify_output_files(clair3, 'VCF')

    def test_clair3_ont(self) -> None:
        """
        actually testing BTyper
        """
        clair3 = Clair3(self.camel)
        clair3.add_input_files({'FASTA': [TestClair3.FILE_FASTA_ONT], 'BAM': [TestClair3.FILE_BAM_ONT]})
        clair3.update_parameters(output_path=self.running_dir, platform='ont',
                                 model_path='/usr/local/bin/lmod/clair3/0.1.12/bin/models/ont/')
        clair3.run(self.running_dir)
        self.verify_output_files(clair3, 'VCF')


if __name__ == '__main__':
    unittest.main()
