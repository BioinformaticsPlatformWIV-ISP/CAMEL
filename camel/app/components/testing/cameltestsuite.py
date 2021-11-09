import logging
import os
import unittest
from pathlib import Path

import shutil
import tempfile

from camel.app.camel import Camel
from camel.app.tools.tool import Tool


class CamelTestSuite(unittest.TestCase):
    """
    Base class for Camel test suites.
    """
    running_dir = None
    camel = Camel.get_instance()

    @staticmethod
    def get_test_file_dir(*args) -> Path:
        """
        Retrieves thee directory with test files.
        :param args: Sub-directory or directories
        :return: Test files directory
        """
        dir_test = Path(CamelTestSuite.camel.config['testing']['testfiles_dir'], *args)
        if not dir_test.exists() or not dir_test.is_dir():
            raise FileNotFoundError(f"Cannot find test file directory: {dir_test}")
        return dir_test

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = Path(tempfile.mkdtemp(None, 'camel_', CamelTestSuite.camel.config['temp_dir']))
        logging.debug(f"Directory for testing: {self.running_dir}")

    def tearDown(self) -> None:
        """
        Removes the working directory after executing the tests.
        :return: None
        """
        if os.environ.get('CAMEL_KEEP_TEST_DIRS') == '1':
            logging.debug("Keeping working directory (CAMEL_KEEP_TEST_DIRS)")
            return
        if Path(self.running_dir).exists():
            logging.debug(f"Removing working directory: {self.running_dir}")
            shutil.rmtree(self.running_dir)

    def verify_output_files(self, tool: Tool, key: str, nb_files: int = 1) -> None:
        """
        Verifies if the output files with the given key are created correctly.
        :param tool: Tool
        :param key: Output key
        :param nb_files: Number of output files that should be generated
        :return: None
        """
        self.assertIn(key, tool.tool_outputs, f"Key '{key}' missing from tool outputs")
        self.assertEqual(
            nb_files, len(tool.tool_outputs[key]),
            f"Unexpected number of tools outputs found: {len(tool.tool_outputs[key])}, expected {nb_files}")
        for i in range(0, nb_files):
            output_file_path = tool.tool_outputs[key][0].path
            self.assertTrue(output_file_path.exists(), f"Output file '{key}' (index: {i}) does not exist")
            self.assertGreater(output_file_path.stat().st_size, 0, f"Output file '{key}' (index: {i}) is empty")
