import argparse
import unittest

import os
import tempfile

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.scripts.pointfinder.mainpointfinder import MainPointFinder


class TestPointFinder(unittest.TestCase):
    """
    Tests the PointFinder tool.
    """
    camel = Camel()
    running_dir = None

    # Input files
    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'])
    input_fasta_file = ToolIOFile(os.path.join(test_file_dir, 'pointfinder', 'ref_ecoli.fasta'))

    def setUp(self) -> None:
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(prefix='camel_', dir=TestPointFinder.camel.config['temp_dir'])

    def test_pointfinder(self) -> None:
        """
        Tests the PointFinder main script
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            sample_name='test_sample',
            fasta=self.input_fasta_file.path,
            fasta_name=os.path.basename(self.input_fasta_file.basename),
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            working_dir=self.running_dir
        )
        main = MainPointFinder(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)
