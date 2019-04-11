import argparse
import unittest

import os
import tempfile

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.scripts.genedetection.mainmakegenedetectiondb import MainMakeGeneDetectionDB


class TestMakeGeneDetectionDB(unittest.TestCase):
    """
    Tests the gene detection create DB tool.
    """

    camel = Camel()
    running_dir = None

    # Input files
    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'], 'gene_detection')
    input_fasta = ToolIOFile(os.path.join(test_file_dir, 'test_input_db.fasta'))

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(None, 'camel_', TestMakeGeneDetectionDB.camel.config['temp_dir'])

    def test_gene_detection_create_db(self) -> None:
        """
        Tests the gene detection create db main script.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            fasta=TestMakeGeneDetectionDB.input_fasta.path,
            fasta_name=None,
            working_dir=self.running_dir,
            identity_cutoff=85
        )
        main = MainMakeGeneDetectionDB(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_gene_detection_create_db_spaces_in_name(self) -> None:
        """
        Tests the gene detection create db main script.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            fasta=TestMakeGeneDetectionDB.input_fasta.path,
            fasta_name='this is my database.fasta',
            working_dir=self.running_dir,
            identity_cutoff=85
        )
        main = MainMakeGeneDetectionDB(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)


if __name__ == '__main__':
    unittest.main()
