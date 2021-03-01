import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.scripts.pointfinder.mainpointfinder import MainPointFinder
from camel.tests import longRunningTest


class TestPointFinder(CamelTestSuite):
    """
    Tests the PointFinder tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('pointfinder')
    input_fasta_file = test_file_dir / 'ref_ecoli.fasta'
    input_fasta_file_salmonella = test_file_dir / 'salmonella_lt2_ref.fasta'
    input_fastq = [
        test_file_dir / 'reads_illumina_1.fastq',
        test_file_dir / 'reads_illumina_2.fastq',
    ]

    def test_pointfinder_main(self) -> None:
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

    def test_pointfinder_main_with_pubmed_link(self) -> None:
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
    def test_pointfinder_main_fastq_input(self) -> None:
        """
        Tests the PointFinder main script with FASTQ input.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(self.input_fastq[0]), str(self.input_fastq[1]),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--species', 'escherichia_coli',
            '--assembly-cov-cutoff', '2',
            '--assembly-min-contig-length', '500'
        ]
        main = MainPointFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
