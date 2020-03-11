import unittest
from pathlib import Path

import tempfile
import yaml

from camel.app.camel import Camel
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.scripts.stecpipeline import CONFIG_DATA
from camel.scripts.stecpipeline.mainstecpipeline import MainSTECPipeline
from camel.tests import longRunningTest


class TestSTECPipeline(unittest.TestCase):
    """
    Tests for the STEC pipeline.
    """

    camel = Camel()
    running_dir = None

    # Input files
    test_file_dir = Path(camel.config['testing']['testfiles_dir'])
    input_fastq_pe = [
        test_file_dir / 'pipelines' / 'STEC-591_S13-ds_1.fastq.gz',
        test_file_dir / 'pipelines' / 'STEC-591_S13-ds_2.fastq.gz'
    ]

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = Path(tempfile.mkdtemp(None, 'camel_', TestSTECPipeline.camel.config['temp_dir']))

    def test_stec_pipeline_typing_db(self) -> None:
        """
        Checks if the databases for the sequence typing are available.
        :return: None
        """
        from camel.app.tools.pipelines.sequence_typing.locussetmanager import LocusSetManager
        with open(CONFIG_DATA) as handle_in:
            config_data = yaml.safe_load(handle_in)

        for key, path in config_data['sequence_typing'].items():
            # Check if scheme exists
            self.assertGreater(Path(path).stat().st_size, 0)

            # Check if metadata can be loaded
            manager = LocusSetManager(Camel.get_instance())
            manager.add_input_files({'DIR': [ToolIODirectory(path)]})
            manager.run(str(self.running_dir))
            self.assertGreater(len(manager.informs), 0)

    def test_stec_pipeline_gene_detection_db(self):
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
    def test_stec_pipeline(self) -> None:
        """
        Tests the STEC pipeline with all assays except for cgMLST.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        path_summary_out = self.running_dir / 'out' / 'summary.tsv'
        args = [
            '--fastq-pe', str(TestSTECPipeline.input_fastq_pe[0]), str(TestSTECPipeline.input_fastq_pe[1]),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--output-tsv', str(path_summary_out),
            '--working-dir', str(self.running_dir)
        ] + [f"--{a.replace('_', '-')}" for a in MainSTECPipeline.CUSTOM_ANALYSES if a != 'cgmlst']
        main = MainSTECPipeline(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
