import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.seqtk.seqtkconvert import SeqtkConvert

class TestSeqtk(CamelTestSuite):
    """
    Tests for the seqtk module.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('seqtk')
    FILE_FASTA = ToolIOFile(test_file_dir / 'lm1_1.fasta')
    FILE_FASTQ = ToolIOFile(test_file_dir / 'lm1_1.fastq')

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
