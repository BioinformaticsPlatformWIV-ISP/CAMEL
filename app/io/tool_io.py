import abc
from abc import ABCMeta


class ToolIO(object):
    """
    Class that represents the input or output of a tool.
    """

    __metaclass__ = ABCMeta

    def __init__(self):
        """
        Initializes a tool input / output.
        """
        pass

    @abc.abstractmethod
    def is_valid(self):
        """
        Checks whether the tool input / output is valid.
        :return: None
        """
        pass
