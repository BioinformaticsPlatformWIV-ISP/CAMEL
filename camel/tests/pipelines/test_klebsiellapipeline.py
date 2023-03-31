import tempfile
import unittest
from pathlib import Path

import yaml

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.tools.pipelines.sequence_typing.locussetmanager import LocusSetManager
from camel.scripts.klebsiellapipeline import CONFIG_DATA
from camel.scripts.klebsiellapipeline.mainklebsiellapipeline import MainKlebsiellaPipeline
from camel.tests import longRunningTest


class TestKlebsiellaPipeline(unittest.TestCase):
    """
    Tests for the Klebsiella pipeline.
    """

    camel = Camel()
    running_dir = None

    # Input files
    test_file_dir = Path(camel.config['testing']['testfiles_dir'])
    input_fastq_pe = [
        test_file_dir / 'pipelines' / 'Kpneumoniae-SRR4046826-ds_1.fastq.gz',
        test_file_dir / 'pipelines' / 'Kpneumoniae-SRR4046826-ds_2.fastq.gz'
    ]

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = Path(tempfile.mkdtemp(None, 'camel_', TestKlebsiellaPipeline.camel.config['temp_dir']))

    def test_klebsiella_pipeline_typing_db(self) -> None:
        """
        Checks if the databases for the sequence typing are available.
        :return: None
        """
        with open(CONFIG_DATA) as handle_in:
            config_data = yaml.safe_load(handle_in)

        for key, scheme_data in config_data['sequence_typing'].items():
            # Check if scheme exists
            self.assertGreater(Path(scheme_data['path']).stat().st_size, 0)

            # Check if metadata can be loaded
            manager = LocusSetManager(Camel.get_instance())
            manager.add_input_files({'DIR': [ToolIODirectory(Path(scheme_data['path']))]})
            manager.run(self.running_dir)
            self.assertGreater(len(manager.informs), 0)

    @longRunningTest()
    def test_klebsiella_pipeline(self):
        """
        Tests the Klebsiella pipeline.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestKlebsiellaPipeline.input_fastq_pe[0]), str(TestKlebsiellaPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir),
            '--threads', '8'
        ] + [
            f"--{a.replace('_', '-')}" for a in MainKlebsiellaPipeline.CUSTOM_ANALYSES if a not in (
                'cgmlst', 'scgmlst')]
        main = MainKlebsiellaPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
