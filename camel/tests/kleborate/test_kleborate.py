import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.kleborate.kleborate import Kleborate


class TestKleborate(CamelTestSuite):
    """
    Tests the Kleborate tool.
    """
    # Get test file and reference file directories
    test_file_dir = CamelTestSuite.get_test_file_dir('kleborate')
    input_fasta = test_file_dir / 'GCF_002248955.1.fna.gz'

    def test_kleborate(self) -> None:
        """
        Tests the Kleborate tool.
        :return: None
        """
        kleborate = Kleborate(self.camel)
        kleborate.add_input_files({'FASTA': [ToolIOFile(TestKleborate.input_fasta)]})
        kleborate.run(self.running_dir)
        self.verify_output_files(kleborate, 'TXT')

    def test_kleborate_different_confidence(self) -> None:
        """
        Tests the Kleborate tool.
        :return: None
        """
        kleborate = Kleborate(self.camel)
        kleborate.add_input_files({'FASTA': [ToolIOFile(TestKleborate.input_fasta)]})
        kleborate.update_parameters(min_kaptive_confidence='Low')
        kleborate.run(self.running_dir)
        self.verify_output_files(kleborate, 'TXT')

    def test_kleborate_resistance(self) -> None:
        """
        Tests the Kleborate tool.
        :return: None
        """
        kleborate = Kleborate(self.camel)
        kleborate.add_input_files({'FASTA': [ToolIOFile(TestKleborate.input_fasta)]})
        kleborate.update_parameters(resistance=True)
        kleborate.run(self.running_dir)
        self.verify_output_files(kleborate, 'TXT')

    def test_kleborate_all(self) -> None:
        """
        Tests the Kleborate tool.
        :return: None
        """
        kleborate = Kleborate(self.camel)
        kleborate.add_input_files({'FASTA': [ToolIOFile(TestKleborate.input_fasta)]})
        kleborate.update_parameters(all=True)
        kleborate.run(self.running_dir)
        self.verify_output_files(kleborate, 'TXT')


if __name__ == '__main__':
    unittest.main()
