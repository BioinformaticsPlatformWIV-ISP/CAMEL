import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.resfinder.resfinder import ResFinder
from camel.scripts.resfinder.mainresfinder import MainResFinder


class TestResFinder(CamelTestSuite):
    """
    Initializes this testing tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('resfinder')
    FILE_FASTA_1 = ToolIOFile(test_file_dir / 'ref_ecoli.fasta')
    FILE_FASTA_2 = ToolIOFile(test_file_dir / 'salmonella_lt2_ref.fasta')
    FILE_FASTA_3 = ToolIOFile(test_file_dir / 'assembly-VAR305.fasta')
    FILE_FASTQ_1 = ToolIOFile(test_file_dir / 'reads_illumina_1.fastq')
    FILE_FASTQ_2 = ToolIOFile(test_file_dir / 'reads_illumina_2.fastq')
    FILE_FASTQ_ENTERO_1 = ToolIOFile(test_file_dir / 'reads_entero_1.fastq.gz')
    FILE_FASTQ_ENTERO_2 = ToolIOFile(test_file_dir / 'reads_entero_2.fastq.gz')
    DB_RESFINDER = Path(CamelTestSuite.camel.config['db_root'], 'resfinder4')

    def test_resfinder_main_fasta(self) -> None:
        """
        Tests the ResFinder main script.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(self.FILE_FASTA_3),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--db-directory', str(self.DB_RESFINDER),
            '--acquired',
            '--acq-overlap', '40',
            '--point',
            '--min-cov', '60',
            '--threshold', '80',
            '--species', 'Staphylococcus_aureus'
        ]
        main = MainResFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_resfinder_main_fastq(self) -> None:
        """
        Tests the ResFinder main script with FASTQ files.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(self.FILE_FASTQ_1), str(self.FILE_FASTQ_2),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--db-directory', str(self.DB_RESFINDER),
            '--point',
            '--acquired',
            '--acq-overlap', '45',
            '--min-cov', '60',
            '--threshold', '80',
            '--species', 'Escherichia_coli'
        ]
        main = MainResFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_resfinder_fasta(self) -> None:
        """
        Actually testing ResFinder with contigs file.
        :return: None
        """
        resfinder = ResFinder(self.camel)
        resfinder.add_input_files({'FASTA': [TestResFinder.FILE_FASTA_1], 'DIR': [ToolIODirectory(self.DB_RESFINDER)]})
        resfinder.update_parameters(output_path=self.running_dir, min_cov=0.6, threshold=0.8, acquired=True)
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV_genes')
        self.verify_output_files(resfinder, 'TSV_pheno_general')

    def test_resfinder_fastq(self) -> None:
        """
        Testing resfinder with paired-end fastq reads.
        :return: None
        """
        resfinder = ResFinder(self.camel)
        resfinder.add_input_files({'FASTQ_PE': [TestResFinder.FILE_FASTQ_1, TestResFinder.FILE_FASTQ_2],
                                   'DIR': [ToolIODirectory(self.DB_RESFINDER)]})
        resfinder.update_parameters(output_path=self.running_dir, min_cov=0.6, threshold=0.8, acquired=True)
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV_genes')
        self.verify_output_files(resfinder, 'TSV_pheno_general')

    def test_resfinder_pointfinder_fasta(self) -> None:
        """
        Testing resfinder with pointfinder mode and fasta file.
        :return: None
        """
        resfinder = ResFinder(self.camel)
        resfinder.add_input_files({'FASTA': [TestResFinder.FILE_FASTA_1], 'DIR': [ToolIODirectory(self.DB_RESFINDER)]})
        resfinder.update_parameters(output_path=self.running_dir, min_cov=0.6, threshold=0.8, point=True,
                                    species='"escherichia coli"')
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV_point')
        self.verify_output_files(resfinder, 'TSV_pheno_general')
        self.verify_output_files(resfinder, 'TSV_pheno_species')

    def test_resfinder_pointfinder_fastq(self) -> None:
        """
        Testing resfinder with pointfinder mode and fastq files.
        :return: None
        """
        resfinder = ResFinder(self.camel)
        resfinder.add_input_files({'FASTQ_PE': [TestResFinder.FILE_FASTQ_1, TestResFinder.FILE_FASTQ_2],
                                   'DIR': [ToolIODirectory(self.DB_RESFINDER)]})
        resfinder.update_parameters(output_path=self.running_dir, min_cov=0.6, threshold=0.8, point=True,
                                    species='"escherichia coli"')
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV_point')
        self.verify_output_files(resfinder, 'TSV_pheno_general')
        self.verify_output_files(resfinder, 'TSV_pheno_species')

    def test_resfinder_fastq_enterococcus(self) -> None:
        """
        Testing resfinder with enterococcus test fastq files.
        :return: None
        """
        resfinder = ResFinder(self.camel)
        resfinder.add_input_files({'FASTQ_PE': [TestResFinder.FILE_FASTQ_ENTERO_1, TestResFinder.FILE_FASTQ_ENTERO_2]})
        resfinder.update_parameters(output_path=self.running_dir, min_cov=0.9, threshold=0.8, point=True,
                                    species='"enterococcus faecalis"')
        self.verify_output_files(resfinder, 'TSV_point')
        self.verify_output_files(resfinder, 'TSV_pheno_general')
        self.verify_output_files(resfinder, 'TSV_pheno_species')


if __name__ == '__main__':
    unittest.main()
