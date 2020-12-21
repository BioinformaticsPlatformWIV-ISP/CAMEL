import unittest
from pathlib import Path

import yaml

from camel.app.camel import Camel
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.scripts.enterococcuspipeline import CONFIG_DATA
from camel.scripts.enterococcuspipeline.mainenterococcuspipeline import MainEnterococcusPipeline
from camel.tests import longRunningTest


class TestEnterococcusPipeline(CamelTestSuite):
    """
    Tests for the Enterococcus pipeline.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir()
    input_faecalis_fastq_pe = [
        test_file_dir / 'pipelines' / 'Enterococcus_faecalis-SRR12362697-ds_1.fastq.gz',
        test_file_dir / 'pipelines' / 'Enterococcus_faecalis-SRR12362697-ds_2.fastq.gz'
    ]
    input_faecium_fastq_pe = [
        test_file_dir / 'pipelines' / 'Enterococcus_faecium-SRR12388968-ds_1.fastq.gz',
        test_file_dir / 'pipelines' / 'Enterococcus_faecium-SRR12388968-ds_2.fastq.gz'
    ]

    def test_enterococcus_gene_detection_db(self):
        """
        Checks if the databases for the gene detection are available.
        :return: None
        """
        from camel.app.tools.pipelines.genedetection.dbmanager import DBManager
        with open(CONFIG_DATA) as handle_in:
            config_data = yaml.safe_load(handle_in)

        for key, db_data in config_data['gene_detection'].items():
            # Check if scheme exists
            self.assertGreater(Path(db_data['path']).stat().st_size, 0)

            # Check if metadata and FASTA files can be loaded
            manager = DBManager(Camel.get_instance())
            manager.add_input_files({'DIR': [ToolIODirectory(db_data['path'])]})
            manager.run(str(self.running_dir))
            self.assertGreater(len(manager.tool_outputs), 0)
            self.assertGreater(len(manager.informs), 0)

    @longRunningTest()
    def test_enterococcus_pipeline_faecalis_blast(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestEnterococcusPipeline.input_faecalis_fastq_pe[0]),
            str(TestEnterococcusPipeline.input_faecalis_fastq_pe[1]),
            '--species', 'faecalis',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ] + [f"--{a.replace('_', '-')}" for a in MainEnterococcusPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainEnterococcusPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_enterococcus_pipeline_faecalis_srst2(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestEnterococcusPipeline.input_faecalis_fastq_pe[0]),
            str(TestEnterococcusPipeline.input_faecalis_fastq_pe[1]),
            '--species', 'faecalis',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'srst2'
        ] + [f"--{a.replace('_', '-')}" for a in MainEnterococcusPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainEnterococcusPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_enterococcus_pipeline_faecalis_kma(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestEnterococcusPipeline.input_faecalis_fastq_pe[0]),
            str(TestEnterococcusPipeline.input_faecalis_fastq_pe[1]),
            '--species', 'faecalis',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma'
        ] + [f"--{a.replace('_', '-')}" for a in MainEnterococcusPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainEnterococcusPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_enterococcus_pipeline_faecium_blast(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestEnterococcusPipeline.input_faecium_fastq_pe[0]),
            str(TestEnterococcusPipeline.input_faecium_fastq_pe[1]),
            '--species', 'faecium',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ] + [f"--{a.replace('_', '-')}" for a in MainEnterococcusPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainEnterococcusPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_enterococcus_pipeline_faecium_srst2(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestEnterococcusPipeline.input_faecium_fastq_pe[0]),
            str(TestEnterococcusPipeline.input_faecium_fastq_pe[1]),
            '--species', 'faecium',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'srst2'
        ] + [f"--{a.replace('_', '-')}" for a in MainEnterococcusPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainEnterococcusPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_enterococcus_pipeline_faecium_kma(self) -> None:
        """
        Tests the Neisseria pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestEnterococcusPipeline.input_faecium_fastq_pe[0]),
            str(TestEnterococcusPipeline.input_faecium_fastq_pe[1]),
            '--species', 'faecium',
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma'
        ] + [f"--{a.replace('_', '-')}" for a in MainEnterococcusPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainEnterococcusPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
