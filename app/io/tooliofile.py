import os

from app.components.files.fileutils import FileUtils
from app.io.toolio import ToolIO


class ToolIOFile(ToolIO):
    """
    Class that represents an input / output file of a tool.
    """
    TYPE_NAME = 'file'

    def __init__(self, path, logged=True):
        """
        Initializes a tool input / output file.
        :param path: Path to the file
        :param logged: If True, the output can be logged
        """
        super(ToolIOFile, self).__init__(logged)
        self._path = os.path.abspath(path)

    def __str__(self):
        """
        String representation
        :return: String representation
        """
        return self._path

    def __repr__(self):
        """
        Internal representation
        :return: Internal representation
        """
        return 'ToolIOFile("{}")'.format(self.path)

    def is_valid(self):
        """
        Checks if the tool input / output file is valid.
        :return: True if valid
        """
        if not self.exists:
            return False
        if not os.path.isfile(self.path):
            return False
        return True

    @property
    def path(self):
        """
        Returns the path to the input / output file.
        :return: Path
        """
        return self._path

    @property
    def basename(self):
        """
        Returns the basename of the input / output file.
        :return: Basename
        """
        return os.path.basename(self.path)

    @property
    def file_extension(self):
        """
        Returns the file extension.
        :return: File extension
        """
        return os.path.splitext(self.path)[-1]

    @property
    def exists(self):
        """
        Checks whether this file exists.
        :return: True if the file exists, False otherwise
        """
        return os.path.isfile(self._path)

    @property
    def size(self):
        """
        Returns the size of this file.
        :return: Size
        """
        return os.path.getsize(self._path)

    @property
    def hash(self):
        """
        Returns the hash value of this file.
        :return: Hash
        """
        return FileUtils.hash_file(self.path)
