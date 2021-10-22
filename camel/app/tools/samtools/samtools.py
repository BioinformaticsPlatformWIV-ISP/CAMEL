import abc

from camel.app.camel import Camel
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.tools.toolpipeable import ToolPipeable


class Samtools(ToolPipeable, metaclass=abc.ABCMeta):
    """
    Super class for samtools.
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

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        super(Samtools, self)._execute_tool()

    def _check_stderr(self) -> None:
        """
        Validate if the program ran correctly by checking the standard error.
        :return: None
        """
        if any(keyword in self._command.stderr.lower() for keyword in ('aborted', 'error')):
            raise ToolExecutionError(f"{self.name} failed: '{self._command.stderr}'")
