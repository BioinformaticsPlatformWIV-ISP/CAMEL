import abc

from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.tools.tool import Tool


class Filter(Tool, metaclass=abc.ABCMeta):
    """
    Base class for variant filters.
    """

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        if 'VCF_GZ' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No compressed VCF input found")
        super(Filter, self)._check_input()
