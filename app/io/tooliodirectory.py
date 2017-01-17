import os

from app.components.files.fileutils import FileUtils
from app.io.toolio import ToolIO


class ToolIODirectory(ToolIO):
    """
    Class that represents an input / output directory of a tool.
    """
    TYPE_NAME = 'dir'

    def __init__(self, path, logged=True):
        """
        Initializes a tool input / output directory.
        :param path: Path to the directory
        :param logged: If True, the output can be logged
        """
        super(ToolIODirectory, self).__init__(logged)
        self._path = path

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
        return 'ToolIODirectory("{}")'.format(self.path)

    def is_valid(self):
        """
        Checks if the tool input / output directory is valid.
        :return: True if valid
        """
        if not self.exists:
            return False
        if not os.path.isdir(self.path):
            return False
        return True

    @property
    def hash(self):
        """
        Returns the hash value.
        :return: Hash value
        """
        return FileUtils.hash_directory(self.path)

    @property
    def path(self):
        """
        Returns the path to the input / output directory.
        :return: Path
        """
        return self._path

    @property
    def basename(self):
        """
        Returns the basename of the input / output directory.
        :return: Basename
        """
        return os.path.basename(self.path)

    @property
    def exists(self):
        """
        Checks whether this directory exists.
        :return: True if the directory exists, False otherwise
        """
        return os.path.isdir(self._path)
