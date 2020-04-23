import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.scripts.sccmecfinder.mainsccmecfinder import MainSCCmecFinder


class TestSCCmecFinder(CamelTestSuite):
    """
    Tests the SCCmecFinder tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('sccmec')
    input_fasta_file = test_file_dir / 'MSSA476.fasta'

    def test_sccmecfinder(self) -> None:
        """
        Tests the SCCmecFinder tool
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(TestSCCmecFinder.input_fasta_file),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--db-mec-genes', '/db/gene_detection/SCCmec_genes',
            '--profiles-mec-genes', '/db/pipelines/saureus/SCCmec_genes/profiles.yml'
        ]
        main = MainSCCmecFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
