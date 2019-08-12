import unittest
from pathlib import Path

import os
import tempfile

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.minimap2.minimap2mapping import Minimap2Mapping


class TestMinimap2(unittest.TestCase):
    """
    Tests the minimap2 tool.
    """

    camel = Camel()
    running_dir = None

    test_file_dir = Path(camel.config['testing']['testfiles_dir']) / 'minion'
    FILE_REF_GENOME = ToolIOFile(test_file_dir / 'NC_002695.2.fasta')
    FILE_FASTQ = ToolIOFile(test_file_dir / 'fastq_minion_stec.fastq')

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(None, 'camel_', TestMinimap2.camel.config['temp_dir'])

    def test_minimap2(self):
        """
        Tests the minimap2 tool.
        :return: None
        """
        minimap2 = Minimap2Mapping(self.camel)
        minimap2.add_input_files({
            'FASTA': [TestMinimap2.FILE_REF_GENOME],
            'FASTQ': [TestMinimap2.FILE_FASTQ]
        })
        minimap2.run(self.running_dir)
        self.assertTrue('SAM' in minimap2.tool_outputs, "No VCF output generated")
        output_file = minimap2.tool_outputs['SAM'][0].path
        self.assertTrue(os.path.isfile(output_file))
        self.assertGreater(os.path.getsize(output_file), 0)


if __name__ == '__main__':
    unittest.main()
