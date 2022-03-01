import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.fqtools.fqstats import FqStats


class TestFqtools(CamelTestSuite):
    """
    Tests the fqtools tool suite.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('fqstats')
    FILE_FQ_FWD = ToolIOFile(test_file_dir / 'ecoli_1.fastq')
    FILE_FQ_REV = ToolIOFile(test_file_dir / 'ecoli_2.fastq')
    FILE_FQ_GZ_FWD = ToolIOFile(test_file_dir / 'ecoli_1.fastq.gz')
    FILE_FQ_GZ_REV = ToolIOFile(test_file_dir / 'ecoli_2.fastq.gz')

    def test_fqtools_stats(self) -> None:
        """
        Tests the fqtools stats function.
        :return: None
        """
        fqtools_stats = FqStats(self.camel)
        fqtools_stats.add_input_files({
            'FASTQ': [TestFqtools.FILE_FQ_FWD, TestFqtools.FILE_FQ_REV]})
        fqtools_stats.run(self.running_dir)
        self.assertIn('stats', fqtools_stats.informs)
        self.assertEquals(len(fqtools_stats.informs['stats']), 2)

    def test_fqtools_stats_gzip_input(self) -> None:
        """
        Tests the fqtools stats function with gzipped input.
        :return: None
        """
        fqtools_stats = FqStats(self.camel)
        fqtools_stats.add_input_files({
            'FASTQ': [TestFqtools.FILE_FQ_GZ_FWD, TestFqtools.FILE_FQ_GZ_REV]})
        fqtools_stats.run(self.running_dir)
        self.assertIn('stats', fqtools_stats.informs)
        self.assertEquals(len(fqtools_stats.informs['stats']), 2)


if __name__ == '__main__':
    unittest.main()
