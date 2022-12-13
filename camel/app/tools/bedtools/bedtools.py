import abc

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.tools.tool import Tool


class Bedtools(Tool, metaclass=abc.ABCMeta):

    """
    The master class for Bedtools toolset
    """

    def __init__(self, tool_name: str, version: str, camel: Camel) -> None:
        """
        Initialize a samtools tool.
        :param tool_name: Tool name
        :param version: Tool version
        :param camel: Camel instance
        :return: None
        """
        super().__init__(tool_name, version,  camel)
        self._required_inputs = []

    def _check_command_output(self) -> None:
        """
        Validate if the program ran correctly by checking the standard error.
        :return: None
        """
        if any(keyword in self.stderr.lower() for keyword in ('aborted', 'error')):
            raise ToolExecutionError(f'{self.name} failed: {self.stderr}')

    def _check_required_inputs(self) -> None:
        """
        Check required input.
        :return: None
        """
        if len(self._required_inputs) == 0:
            return
        try:
            next(key for key in self._required_inputs if key in self._tool_inputs)
        except StopIteration:
            raise InvalidInputSpecificationError(f"Required inputs missing ({', '.join(self._required_inputs)})")
