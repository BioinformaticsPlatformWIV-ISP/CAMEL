import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.seqtk.seqtksize import SeqtkSize


class TestSeqtk(CamelTestSuite):
    """
    Initializes the seqtk tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('seqtk')
    FILE_FASTQ_A = ToolIOFile(test_file_dir / 'lm1_1.fastq')
    FILE_FASTQ_B = ToolIOFile(test_file_dir / 'lm1_2.fastq')

    def test_seqtk_size(self) -> None:
        """
        Testing SeqKit seq with contigs file.
        """
        seqtk_size = SeqtkSize(self.camel)
        seqtk_size.add_input_files({'FASTQ': [TestSeqtk.FILE_FASTQ_A, TestSeqtk.FILE_FASTQ_B]})
        seqtk_size.run(self.running_dir)
        self.assertIn('stats', seqtk_size.informs)
        self.assertEqual(len(seqtk_size.informs['stats']), 2)


if __name__ == '__main__':
    unittest.main()
