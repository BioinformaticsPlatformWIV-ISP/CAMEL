import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.scripts.pointfinder.mainpointfinder import MainPointFinder
from camel.tests import longRunningTest


class TestPointFinder(CamelTestSuite):
    """
    Tests the PointFinder tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir()
    input_fasta_file = test_file_dir / 'pointfinder' / 'ref_ecoli.fasta'
    input_fasta_file_salmonella = test_file_dir / 'pointfinder' / 'salmonella_lt2_ref.fasta'
    input_fastq_raw_galaxy = [
        test_file_dir / 'workflows' / 'dataset_fwd_11.dat',
        test_file_dir / 'workflows' / 'dataset_rev_10.dat'
    ]

    def test_pointfinder(self) -> None:
        """
        Tests the PointFinder main script
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta_file),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--species', 'escherichia_coli'
        ]
        main = MainPointFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_pointfinder_with_pubmed_link(self) -> None:
        """
        Tests the PointFinder main script with a link to PubMed.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta_file_salmonella),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--species', 'salmonella'
        ]
        main = MainPointFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    @longRunningTest()
    def test_pointfinder_fastq_input(self) -> None:
        """
        Tests the PointFinder main script with FASTQ input.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(self.input_fastq_raw_galaxy[0]), str(self.input_fastq_raw_galaxy[1]),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--species', 'escherichia_coli',
            '--assembly-cov-cutoff', '5',
            '--assembly-min-contig-length', '500'
        ]
        main = MainPointFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
