from app.components.files.fileutils import FileUtils
from app.io.toolio import ToolIO


class ToolIOValue(ToolIO):
    """
    Class that represents an input / output value of a tool.
    """
    TYPE_NAME = 'value'

    def __init__(self, value, logged=True):
        """
        Initializes a tool input / output value.
        :param value: Value
        :param logged: If True, the output can be logged
        """
        super(ToolIOValue, self).__init__(logged)
        self._value = value

    @property
    def value(self):
        """
        Returns the value.
        :return: Value
        """
        return self._value

    @property
    def hash(self):
        """
        Returns the hash value.
        :return: Hash value
        """
        return FileUtils.hash_value(self.value)

    def __str__(self):
        """
        String representation
        :return: String representation
        """
        return str(self.value)

    def __repr__(self):
        """
        Internal representation
        :return: Internal representation representation
        """
        return 'ToolIOValue({})'.format(repr(self.value))

    def is_valid(self):
        """
        Checks if the tool input / output value is valid.
        :return: True if valid
        """
        if self._value is None:
            return False
        return True
