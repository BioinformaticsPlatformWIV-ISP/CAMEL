import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.resfinder.resfinder import ResFinder
from camel.scripts.resfinder.mainresfinder import MainResFinder


class TestResFinder(CamelTestSuite):
    """
    Initializes this testing tool
    """

    test_file_dir = Path('/testdata/camel/pointfinder/')
    FILE_FASTA_1 = ToolIOFile(test_file_dir / 'ref_ecoli.fasta')
    FILE_FASTA_2 = ToolIOFile(test_file_dir / 'salmonella_lt2_ref.fasta')
    FILE_FASTQ_1 = ToolIOFile(test_file_dir / 'reads_illumina_1.fastq')
    FILE_FASTQ_2 = ToolIOFile(test_file_dir / 'reads_illumina_2.fastq')

    def test_resfinder_main(self) -> None:
        """
        Tests the ResFinder main script
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(self.FILE_FASTA_1),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--acquired',
            '--min_cov', '0.6',
            '--threshold', '0.8'
        ]
        main = MainResFinder(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_resfinder_fasta(self) -> None:
        """
        actually testing ResFinder with contigs file
        """
        resfinder = ResFinder(self.camel)
        resfinder.add_input_files({'FASTA': [TestResFinder.FILE_FASTA_1]})
        resfinder.update_parameters(output_path=self.running_dir, min_cov=0.6, threshold=0.8)
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV')

    def test_resfinder_fastq(self) -> None:
        """
        testing resfinder with paired-end fastq reads
        """
        resfinder = ResFinder(self.camel)
        resfinder.add_input_files({'FASTQ_PE': [TestResFinder.FILE_FASTQ_1, TestResFinder.FILE_FASTQ_2]})
        resfinder.update_parameters(output_path=self.running_dir, min_cov=0.6, threshold=0.8)
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV')

    def test_resfinder_pointfinder_fasta(self) -> None:
        """
        testing resfinder with pointfinder mode and fasta file
        """
        resfinder = ResFinder(self.camel)
        resfinder.add_input_files({'FASTA': [TestResFinder.FILE_FASTA_1]})
        resfinder.update_parameters(output_path=self.running_dir, min_cov=0.6, threshold=0.8, point=True,
                                    species='ecoli')
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV')

    def test_resfinder_pointfinder_fastq(self) -> None:
        """
        testing resfinder with pointfinder mode and fastq files
        """
        resfinder = ResFinder(self.camel)
        resfinder.add_input_files({'FASTQ_PE': [TestResFinder.FILE_FASTQ_1, TestResFinder.FILE_FASTQ_2]})
        resfinder.update_parameters(output_path=self.running_dir, min_cov=0.6, threshold=0.8, point=True,
                                    species='ecoli')
        resfinder.run(self.running_dir)
        self.verify_output_files(resfinder, 'TSV')


if __name__ == '__main__':
    unittest.main()
