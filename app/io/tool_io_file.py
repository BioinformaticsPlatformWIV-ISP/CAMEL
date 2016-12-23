import os

from app.io.tool_io import ToolIO


class ToolIOFile(ToolIO):
    """
    Class that represents an input / output file of a tool.
    """

    def __init__(self, path):
        """
        Initializes a tool input / output file.
        :param path: Path to the file
        """
        super(ToolIO, self).__init__()
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
        :return: Internal representation representation
        """
        return 'ToolIOFile("{}")'.format(self.path)

    def is_valid(self):
        """
        Checks if the tool input / output file is valid.
        :return: True if valid
        """
        if not self.exists:
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
