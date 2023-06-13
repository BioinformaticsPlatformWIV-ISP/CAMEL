import unittest
from pathlib import Path

from camel.app.camel import Camel
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.tools.pipelines.sequence_typing.locussetmanager import LocusSetManager
from camel.scripts.bacilluspipeline.mainbacilluspipeline import MainBacillusPipeline
from camel.tests import longRunningTest


class TestBacillusPipeline(CamelTestSuite):
    """
    Tests for the Bacillus pipeline.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir()
    input_fastq_pe_subtilis = [
        test_file_dir / 'pipelines' / 'Bsubtilis-SRR10568181_1.fastq.gz',
        test_file_dir / 'pipelines' / 'Bsubtilis-SRR10568181_2.fastq.gz'
    ]
    input_fastq_pe_cereus = [
        test_file_dir / 'pipelines' / 'Bcereus-SRR7067969-ds_1.fastq.gz',
        test_file_dir / 'pipelines' / 'Bcereus-SRR7067969-ds_2.fastq.gz'
    ]

    def test_bacillus_pipeline_typing_db(self) -> None:
        """
        Checks if the databases for the sequence typing are available.
        :return: None
        """
        sequence_typing_dict = {'mlst_cereus': {'path': '/db/sequence_typing/bacillus_cereus/mlst'},
                                'mlst_subtilis': {'path': '/db/sequence_typing/bacillus_subtilis/mlst'}}
        for key, scheme_data in sequence_typing_dict.items():
            # Check if scheme exists
            self.assertGreater(Path(scheme_data['path']).stat().st_size, 0)

            # Check if metadata can be loaded
            manager = LocusSetManager(Camel.get_instance())
            manager.add_input_files({'DIR': [ToolIODirectory(Path(scheme_data['path']))]})
            manager.run(self.running_dir)
            self.assertGreater(len(manager.informs), 0)

    @longRunningTest()
    def test_bacillus_subtilis_pipeline_blast(self) -> None:
        """
        Tests the Bacillus pipeline with blast based detection.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'
        args = [
                   '--fastq-pe', str(TestBacillusPipeline.input_fastq_pe_subtilis[0]),
                   str(TestBacillusPipeline.input_fastq_pe_subtilis[1]),
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir),
                   '--detection-method', 'blast',
                   '--threads', '8',
                   '--species', 'subtilis'
               ] + [f"--{a.replace('_', '-')}" for a in MainBacillusPipeline.CUSTOM_ANALYSES if 'cgmlst' not in a]
        main = MainBacillusPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_bacillus_cereus_pipeline_blast(self) -> None:
        """
        Tests the Bacillus pipeline with blast based detection.
        :return: None
        """
        path_report_out = Path(self.running_dir) / 'out' / 'report.html'
        path_summary_out = Path(self.running_dir) / 'out' / 'summary.tsv'
        args = [
                   '--fastq-pe', str(TestBacillusPipeline.input_fastq_pe_cereus[0]),
                   str(TestBacillusPipeline.input_fastq_pe_cereus[1]),
                   '--output-html', str(path_report_out),
                   '--output-dir', str(path_report_out.parent),
                   '--output-tsv', str(path_summary_out),
                   '--working-dir', str(self.running_dir),
                   '--detection-method', 'blast',
                   '--threads', '8',
                   '--species', 'cereus'
               ] + [f"--{a.replace('_', '-')}" for a in MainBacillusPipeline.CUSTOM_ANALYSES if 'cgmlst' not in a]
        main = MainBacillusPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
