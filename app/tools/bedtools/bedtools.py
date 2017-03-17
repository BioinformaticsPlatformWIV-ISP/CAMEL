import abc

from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.error.toolexecutionerror import ToolExecutionError
from app.tools.tool import Tool


class Bedtools(Tool):

    """
    The master class for Bedtools toolset
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, tool_name, version, camel):
        """
        Initialize a samtools tool.
        :param tool_name: Tool name
        :param version: Tool version
        :param camel: Camel instance
        :return: None
        """
        super(Bedtools, self).__init__(tool_name, version,  camel)
        self._required_inputs = []

    def _check_command_output(self):
        """
        Validate if the program ran correctly by checking the standard error.
        :return: None
        """
        if any(keyword in self.stderr.lower() for keyword in ('aborted', 'error')):
            raise ToolExecutionError("{!r} failed, stderr: {}".format(self.name, self.stderr))

    def _check_required_inputs(self):
        """
        Check required input
        :return: None
        """
        if self._required_inputs:
            for input_type in self._required_inputs:
                if input_type not in self._tool_inputs:
                    raise InvalidInputSpecificationError("No required {!r} input found".format(input_type))
