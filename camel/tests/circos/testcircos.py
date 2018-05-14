import unittest

import os
import tempfile

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.circos.circos import Circos


class TestCircos(unittest.TestCase):
    """
    Tests Circos.
    """

    camel = Camel()
    running_dir = None
    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'], 'circos')

    FILE_CIRCOS_CONFIG = ToolIOFile(os.path.join(test_file_dir, 'hello_world.txt'))

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(None, 'camel_', self.camel.config['temp_dir'])

    def test_circos(self):
        """
        Tests circos.
        :return: None
        """
        circos = Circos(self.camel)
        circos.add_input_files({'TXT': [TestCircos.FILE_CIRCOS_CONFIG]})
        circos.run(self.running_dir)
        self.assertTrue(len(circos.tool_outputs) > 0, "No outputs generated")


if __name__ == '__main__':
    unittest.main()
