import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.scripts.mockpipeline.mainmockpipeline import MainMockPipeline


class TestMocksPipeline(CamelTestSuite):
    """
    Tests for the mock pipeline.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('hybrid_assembly')

    def test_mock_pipeline_illumina(self) -> None:
        """
        Tests the mock pipeline with Illumina input data.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'

        pipeline = MainMockPipeline([
            '--fastq-pe',
            str(TestMocksPipeline.test_file_dir / 'ncbi_region_1.fastq'),
            str(TestMocksPipeline.test_file_dir / 'ncbi_region_2.fastq'),
            '--read-type', 'illumina',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8',
        ])
        pipeline.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_mock_pipeline_ont(self) -> None:
        """
        Tests the mock pipeline with ONT input data.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'

        pipeline = MainMockPipeline([
            '--fastq-se',
            str(TestMocksPipeline.test_file_dir / 'ncbi_region_ont.fastq.gz'),
            '--read-type', 'nanopore',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8',
        ])
        pipeline.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    # def test_mock_pipeline_hybrid(self) -> None:
    #     """
    #     Tests the mock pipeline with hybrid input data (Illumina + ONT).
    #     :return: None
    #     """
    #     path_report_out = Path(self.running_dir) / 'out' / 'report.html'
    #     path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'
    #
    #     pipeline = MainMockPipeline([
    #         '--fastq-pe', str(TestMocksPipeline.input_fastq_se_subtilis),
    #         '--read-type', 'nanopore',
    #         '--output-html', str(path_report_out),
    #         '--output-dir', str(path_report_out.parent),
    #         '--output-tsv', str(path_summary_out),
    #         '--working-dir', str(self.running_dir),
    #         '--detection-method', 'blast',
    #         '--threads', '8',
    #     ])
    #     pipeline.run()
    #     self.assertGreater(path_report_out.stat().st_size, 0)
    #
    # def test_mock_pipeline_fasta(self) -> None:
    #     """
    #     Tests the mock pipeline with FASTA input data.
    #     :return: None
    #     """
    #     path_report_out = Path(self.running_dir) / 'out' / 'report.html'
    #     path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'
    #
    #     pipeline = MainMockPipeline([
    #         '--fastq-pe', str(TestMocksPipeline.input_fastq_se_subtilis),
    #         '--read-type', 'nanopore',
    #         '--output-html', str(path_report_out),
    #         '--output-dir', str(path_report_out.parent),
    #         '--output-tsv', str(path_summary_out),
    #         '--working-dir', str(self.running_dir),
    #         '--detection-method', 'blast',
    #         '--threads', '8',
    #     ])
    #     pipeline.run()
    #     self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
