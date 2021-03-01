import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.scripts.resfinderlocal.mainresfinderlocal import MainResFinderLocal


class TestResFinderLocal(CamelTestSuite):
    """
    Tests the ResFinder local tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir()
    input_fasta = test_file_dir / 'workflows' / 'NC_002695.1.fasta'
    input_fasta_galaxy = test_file_dir / 'workflows' / 'dataset_12.dat'
    input_reads_raw = [
        test_file_dir / 'gene_detection' / 'illumina' / 'reads_illumina_1.fastq',
        test_file_dir / 'gene_detection' / 'illumina' / 'reads_illumina_2.fastq']
    input_db = '/db/gene_detection/ResFinder'

    def test_resfinder_local(self) -> None:
        """
        Tests the ResFinder local main script.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(TestResFinderLocal.input_fasta),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--resfinder-db', str(TestResFinderLocal.input_db)
        ]
        main = MainResFinderLocal(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_resfinder_local_spaces(self) -> None:
        """
        Tests the ResFinder local main script on an input file with spaces.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(TestResFinderLocal.input_fasta),
            '--fasta-name', '"file with spaces.fasta"',
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--resfinder-db', str(TestResFinderLocal.input_db),
            '--min-percent-identity', '85',
            '--min-percent-coverage', '50'
        ]
        main = MainResFinderLocal(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_resfinder_local_galaxy_input(self) -> None:
        """
        Tests the ResFinder local main script on an input file with spaces.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(TestResFinderLocal.input_fasta_galaxy),
            '--fasta-name', 'SPAdes on data 82 and data 81: contigs (fasta)',
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--resfinder-db', str(TestResFinderLocal.input_db)
        ]
        main = MainResFinderLocal(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_resfinder_illumina_pe(self) -> None:
        """
        Tests the ResFinder local main script.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(TestResFinderLocal.input_reads_raw[0]), str(TestResFinderLocal.input_reads_raw[1]),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--resfinder-db', str(TestResFinderLocal.input_db)
        ]
        main = MainResFinderLocal(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
