import abc
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.tools.tool import Tool


class BWA(Tool, metaclass=abc.ABCMeta):
    """Super class for reads mapping using BWA"""

    def _check_command_output(self):
        """
        Parse stderr message of BWA cmd to check whether it runs successfully
        :return: None
        """
        if any(err in self.stderr.lower() for err in ("error", "fail")):
            raise ToolExecutionError('Command failed: {!r}\n stderr: {}'.format(self._command.command, self.stderr))
