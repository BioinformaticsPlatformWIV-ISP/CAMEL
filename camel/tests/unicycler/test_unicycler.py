import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.unicycler.unicycler import Unicycler
from camel.tests import longRunningTest


class TestUnicycler(CamelTestSuite):
    """
    Contains tests for the Unicycler tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('unicycler')
    FILE_FQ = ToolIOFile(test_file_dir / 'subsampled_ont.fastq.gz')
    FILE_FQ_1 = ToolIOFile(test_file_dir / 'pilsner_subsampled_1.fastq.gz')
    FILE_FQ_2 = ToolIOFile(test_file_dir / 'pilsner_subsampled_2.fastq.gz')

    @longRunningTest()
    def test_unicycler(self) -> None:
        """
        Tests Unicycler 0.5.0.
        """
        unicycler = Unicycler(self.camel)
        unicycler.add_input_files({'FASTQ_SE': [TestUnicycler.FILE_FQ],
                                   'FASTQ_PE': [TestUnicycler.FILE_FQ_1, TestUnicycler.FILE_FQ_2]})
        unicycler.run(self.running_dir)
        self.verify_output_files(unicycler, 'FASTA')


if __name__ == '__main__':
    unittest.main()
