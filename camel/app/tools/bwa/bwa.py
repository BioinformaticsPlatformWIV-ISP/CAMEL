import abc

from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.tools.toolpipeable import ToolPipeable


class BWA(ToolPipeable, metaclass=abc.ABCMeta):
    """Super class for reads mapping using BWA"""

    def _check_command_output(self) -> None:
        """
        Parse stderr message of BWA cmd to check whether it runs successfully
        :return: None
        """
        if any(err in self.stderr.lower() for err in ('error', 'fail')):
            raise ToolExecutionError(f'Command failed: {self._command.command}\n stderr: {self.stderr}')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        raise NotImplementedError("Method should be implemented by subclass.")
