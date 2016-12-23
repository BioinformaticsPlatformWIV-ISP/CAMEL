from app.io.tool_io import ToolIO


class ToolIOValue(ToolIO):
    """
    Class that represents an input / output value of a tool.
    """

    def __init__(self, value):
        """
        Initializes a tool input / output value.
        :param value: Value
        """
        super(ToolIO, self).__init__()
        self._value = value

    @property
    def value(self):
        """
        Returns the value.
        :return: Value
        """
        return self._value

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
        return 'ToolIOValue("{}")'.format(self.value)

    def is_valid(self):
        """
        Checks if the tool input / output value is valid.
        :return: True if valid
        """
        if self._value is None:
            return False
        return True
