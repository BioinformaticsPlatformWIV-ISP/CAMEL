import abc

from camel.app.camel import Camel
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.tools.tool import Tool


class Lofreq(Tool, metaclass=abc.ABCMeta):
    """
    Super class for Lofreq.
    """

    def __init__(self, tool_name: str, version: str, camel: Camel) -> None:
        """
        Initialize a Lofreq tool.
        :param tool_name: Tool name
        :param version: Tool version
        :param camel: Camel instance
        :return: None
        """
        super().__init__(tool_name, version,  camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        super()._execute_tool()

    def _check_command_output(self) -> None:
        """
        Validates if the program ran correctly by checking the standard error.
        :return: None
        """
        if 'FATAL' in self._command.stderr:
            raise ToolExecutionError(f"{self.name} failed: '{self._command.stderr}'")
