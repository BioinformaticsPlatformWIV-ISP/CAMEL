import unittest

import os

from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper


class TestFileSystemHelper(unittest.TestCase):
    """
    Tests the Fastq utils module.
    """

    camel = Camel.get_instance()
    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'], 'components')

    def test_is_gzipped_uncompressed(self) -> None:
        """
        Tests the get sample name function for MiSEQ format.
        :return: None
        """
        path = os.path.join(TestFileSystemHelper.test_file_dir, 'fq-file.fastq')
        self.assertFalse(FileSystemHelper.is_gzipped(path))

    def test_is_gzipped_compressed(self) -> None:
        """
        Tests the get sample name function for MiSEQ format.
        :return: None
        """
        path = os.path.join(TestFileSystemHelper.test_file_dir, 'fq-file.fastq.gz')
        self.assertTrue(FileSystemHelper.is_gzipped(path))

    def test_is_gzipped_compressed_no_ext(self) -> None:
        """
        Tests the get sample name function for MiSEQ format.
        :return: None
        """
        path = os.path.join(TestFileSystemHelper.test_file_dir, 'fq-file-no-ext')
        self.assertFalse(FileSystemHelper.is_gzipped(path))

    def test_is_gzipped_uncompressed_no_ext(self) -> None:
        """
        Tests the get sample name function for MiSEQ format.
        :return: None
        """
        path = os.path.join(TestFileSystemHelper.test_file_dir, 'fq-file-no-ext-gzipped')
        self.assertTrue(FileSystemHelper.is_gzipped(path))


if __name__ == '__main__':
    unittest.main()
