import unittest
from pathlib import Path

import yaml

from camel.app.camel import Camel
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.scripts.shigellapipeline import CONFIG_DATA
from camel.scripts.shigellapipeline.mainshigellapipeline import MainShigellaPipeline
from camel.tests import longRunningTest


class TestShigellaPipeline(CamelTestSuite):
    """
    Tests for the Shigella pipeline.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines')
    input_fastq_pe = [
        test_file_dir / 'Shigella-S17BD07654_1.fastq.gz',
        test_file_dir / 'Shigella-S17BD07654_2.fastq.gz'
    ]

    def test_shigella_pipeline_typing_db(self) -> None:
        """
        Checks if the databases for the sequence typing are available.
        :return: None
        """
        from camel.app.tools.pipelines.sequence_typing.locussetmanager import LocusSetManager
        with open(CONFIG_DATA) as handle_in:
            config_data = yaml.safe_load(handle_in)

        for key, scheme_data in config_data['sequence_typing'].items():
            # Check if scheme exists
            self.assertGreater(Path(scheme_data['path']).stat().st_size, 0)

            # Check if metadata can be loaded
            manager = LocusSetManager(Camel.get_instance())
            manager.add_input_files({'DIR': [ToolIODirectory(scheme_data['path'])]})
            manager.run(str(self.running_dir))
            self.assertGreater(len(manager.informs), 0)

    def test_shigella_pipeline_gene_detection_db(self):
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
    def test_shigella_pipeline_blast(self) -> None:
        """
        Tests the Shigella pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestShigellaPipeline.input_fastq_pe[0]), str(TestShigellaPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ] + [f"--{a.replace('_', '-')}" for a in MainShigellaPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainShigellaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    @longRunningTest()
    def test_shigella_pipeline_srst2(self) -> None:
        """
        Tests the Shigella pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestShigellaPipeline.input_fastq_pe[0]), str(TestShigellaPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'srst2'
        ] + [f"--{a.replace('_', '-')}" for a in MainShigellaPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainShigellaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
