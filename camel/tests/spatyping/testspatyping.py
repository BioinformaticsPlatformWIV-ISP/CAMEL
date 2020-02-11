import unittest
from pathlib import Path

import tempfile

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.scripts.spatyping.mainspatyping import MainSpaTyping


class TestSpaTyping(unittest.TestCase):
    """
    Tests the spa typing tool.
    """
    camel = Camel()
    running_dir = None

    # Input files
    test_file_dir = Path(camel.config['testing']['testfiles_dir'])
    input_fasta_file = ToolIOFile(test_file_dir / 'spatyping' / 'saureus_ref.fasta')

    def setUp(self) -> None:
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = Path(tempfile.mkdtemp(prefix='camel_', dir=TestSpaTyping.camel.config['temp_dir']))

    def test_spatyping(self) -> None:
        """
        Tests the spa typing main script
        :return: None
        """
        output_file_report = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(TestSpaTyping.input_fasta_file),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--db-path', '/db/pipelines/saureus'
        ]
        main = MainSpaTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
