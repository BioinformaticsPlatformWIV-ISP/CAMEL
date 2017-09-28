import abc
from abc import ABCMeta


class ToolIO(object, metaclass=ABCMeta):
    """
    Class that represents the input or output of a tool.
    """

    def __init__(self, logged):
        """
        Initializes a tool input / output.
        :param logged: If true, the output can be logged
        """
        self._logged = logged

    @property
    def logged(self):
        """
        Returns True if the output is logged.
        :return: True / False
        """
        return self._logged

    @abc.abstractmethod
    def is_valid(self):
        """
        Checks whether the tool input / output is valid.
        :return: None
        """
        pass

    @property
    @abc.abstractmethod
    def hash(self):
        """
        Returns the hash value.
        :return: Hash value
        """
        pass
