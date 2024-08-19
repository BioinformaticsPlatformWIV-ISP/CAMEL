import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.resfinder.resfinder import ResFinder
from camel.app.tools.resfinder.resfinderreporter import ResFinderReporter
from camel.scripts.resfinder.mainresfinder import MainResFinder


class TestResFinder(CamelTestSuite):
    """
    Initializes this testing tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('resfinder')
    FILE_FASTA_1 = ToolIOFile(test_file_dir / 'ref_ecoli.fasta')
    FILE_FASTA_2 = ToolIOFile(test_file_dir / 'salmonella_lt2_ref.fasta')
    FILE_FASTA_3 = ToolIOFile(test_file_dir / 'assembly-VAR305.fasta')
    FILE_FASTA_4 = ToolIOFile(test_file_dir / 'enterococcus_unkown_amr_muts.fasta')
    FILE_FASTQ_1 = ToolIOFile(test_file_dir / 'reads_illumina_1.fastq')
    FILE_FASTQ_2 = ToolIOFile(test_file_dir / 'reads_illumina_2.fastq')
    FILE_FASTA_ENTERO = ToolIOFile(test_file_dir / 'assembly_entero.fasta')
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

    def test_resfinder_main_fasta_unknown_phenotypes(self) -> None:
        """
        Tests the ResFinder main script with mutations with unknown phenotypes.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(self.FILE_FASTA_4),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--db-directory', str(self.DB_RESFINDER),
            '--acquired',
            '--min-cov', '60',
            '--threshold', '80',
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

    def test_resfinder_main_fastq_kleb(self) -> None:
        """
        Tests the ResFinder main script with FASTQ files for species klebsiella.
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
            '--species', 'Klebsiella'
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
        resfinder.add_input_files({
            'FASTQ_PE': [TestResFinder.FILE_FASTQ_1, TestResFinder.FILE_FASTQ_2],
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
        resfinder.add_input_files({
            'FASTA': [TestResFinder.FILE_FASTA_1],
            'DIR': [ToolIODirectory(self.DB_RESFINDER)]})
        resfinder.update_parameters(
            output_path=self.running_dir, min_cov=0.6, threshold=0.8, point=True, species='"escherichia coli"')
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
        resfinder.add_input_files({
            'FASTQ_PE': [TestResFinder.FILE_FASTQ_1, TestResFinder.FILE_FASTQ_2],
            'DIR': [ToolIODirectory(self.DB_RESFINDER)]})
        resfinder.update_parameters(
            output_path=self.running_dir, min_cov=0.6, threshold=0.8, point=True, species='"escherichia coli"')
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV_point')
        self.verify_output_files(resfinder, 'TSV_pheno_general')
        self.verify_output_files(resfinder, 'TSV_pheno_species')

    def test_resfinder_fasta_enterococcus(self) -> None:
        """
        Testing resfinder with an Enterococcus input FASTA file.
        :return: None
        """
        resfinder = ResFinder(self.camel)
        resfinder.add_input_files({
            'FASTA': [TestResFinder.FILE_FASTA_ENTERO],
            'DIR': [ToolIODirectory(self.DB_RESFINDER)]
        })
        resfinder.update_parameters(
            output_path=self.running_dir, min_cov=0.9, threshold=0.8, point=True, acquired=True,
            species='"enterococcus faecalis"')
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV_point')
        self.verify_output_files(resfinder, 'TSV_genes')
        self.verify_output_files(resfinder, 'TSV_pheno_general')
        self.verify_output_files(resfinder, 'TSV_pheno_species')

        reporter = ResFinderReporter(self.camel)
        reporter.add_input_files({
            'TSV_point': resfinder.tool_outputs['TSV_point'],
            'TSV_pheno_general': resfinder.tool_outputs['TSV_pheno_general'],
            'TSV_pheno_species': resfinder.tool_outputs['TSV_pheno_species'],
            'TSV_genes': resfinder.tool_outputs['TSV_genes']})
        reporter.add_input_informs({'resfinder': resfinder.informs})
        reporter.run(self.running_dir)
        output_section = reporter.tool_outputs['VAL_HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)

    def test_resfinder_fasta_enterococcus_no_species(self) -> None:
        """
        Testing resfinder with an Enterococcus input FASTA file.
        :return: None
        """
        # Run ResFinder
        resfinder = ResFinder(self.camel)
        resfinder.add_input_files({
            'FASTA': [TestResFinder.FILE_FASTA_ENTERO],
            'DIR': [ToolIODirectory(self.DB_RESFINDER)]
        })
        resfinder.update_parameters(
            output_path=self.running_dir, min_cov=0.9, threshold=0.8, point=False, acquired=True)
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV_genes')
        self.verify_output_files(resfinder, 'TSV_pheno_general')

        # Run the reporter
        reporter = ResFinderReporter(self.camel)
        reporter.add_input_files({
            'TSV_pheno_general': resfinder.tool_outputs['TSV_pheno_general'],
            'TSV_genes': resfinder.tool_outputs['TSV_genes']})
        reporter.add_input_informs({'resfinder': resfinder.informs})
        reporter.run(self.running_dir)
        output_section = reporter.tool_outputs['VAL_HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)

    def test_resfinder_main_fasta_enterococcus_no_species(self) -> None:
        """
        Tests the ResFinder main script with the Enterococcus input data.
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(self.FILE_FASTA_ENTERO),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--db-directory', str(self.DB_RESFINDER),
            '--acquired',
            '--acq-overlap', '40',
            '--min-cov', '60',
            '--threshold', '80',
        ]
        main = MainResFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
