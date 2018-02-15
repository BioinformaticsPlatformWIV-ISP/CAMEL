from app.components.files.fileutils import FileUtils
from app.io.toolio import ToolIO
from typing import Hashable


class ToolIOValue(ToolIO):
    """
    Class that represents an input / output value of a tool.
    """

    def __init__(self, value: Hashable, logged: bool=True) -> None:
        """
        Initializes a tool input / output value.
        :param value: Value
        :param logged: If True, the output can be logged
        """
        super(ToolIOValue, self).__init__(logged)
        self._value = value

    @property
    def value(self) -> Hashable:
        """
        Returns the value.
        :return: Value
        """
        return self._value

    @property
    def hash(self) -> str:
        """
        Returns the hash value.
        :return: Hash value
        """
        return FileUtils.hash_value(self.value)

    @property
    def type_name(self) -> str:
        """
        Returns the type of the IO object.
        :return: Type value
        """
        return 'value'

    def __str__(self) -> str:
        """
        String representation
        :return: String representation
        """
        return str(self.value)

    def __repr__(self) -> str:
        """
        Internal representation
        :return: Internal representation representation
        """
        return 'ToolIOValue({})'.format(repr(self.value))

    def is_valid(self) -> bool:
        """
        Checks if the tool input / output value is valid.
        :return: True if valid
        """
        return False if self._value is None else True
