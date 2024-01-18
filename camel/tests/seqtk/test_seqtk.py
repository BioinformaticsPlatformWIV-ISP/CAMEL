import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.seqtk.seqtksize import SeqtkSize
from camel.app.tools.seqtk.seqtkconvert import SeqtkConvert

class TestSeqtk(CamelTestSuite):
    """
    Initializes the seqtk tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('seqtk')
    FILE_FASTQ_A = ToolIOFile(test_file_dir / 'lm1_1.fastq')
    FILE_FASTQ_B = ToolIOFile(test_file_dir / 'lm1_2.fastq')
    FILE_FASTA = ToolIOFile(test_file_dir / 'lm1_1.fasta')

    def test_seqtk_size(self) -> None:
        """
        Testing SeqKit seq with contigs file.
        """
        seqtk_size = SeqtkSize(self.camel)
        seqtk_size.add_input_files({'FASTQ': [TestSeqtk.FILE_FASTQ_A, TestSeqtk.FILE_FASTQ_B]})
        seqtk_size.run(self.running_dir)
        self.assertIn('stats', seqtk_size.informs)
        self.assertEqual(len(seqtk_size.informs['stats']), 2)

    def test_seqtk_convert(self) -> None:
        """
        Testing Seqtk seq -a with fastq file.
        """
        seqtk_convert = SeqtkConvert(self.camel)
        seqtk_convert.add_input_files({'FASTQ': [TestSeqtk.FILE_FASTQ]})
        seqtk_convert.update_parameters(output_file='sequences.fasta')
        seqtk_convert.run(self.running_dir)
        self.verify_output_files(seqtk_convert, 'FASTA')

if __name__ == '__main__':
    unittest.main()
