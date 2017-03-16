import abc

from app.error.toolexecutionerror import ToolExecutionError
from app.tools.tool import Tool


class Bowtie2(Tool):
    """Super class for read mapping using Bowtie2"""
    __metaclass__ = abc.ABCMeta

    def _check_command_output(self):
        """
        Parse stderr message of Bowtie2 cmd to check whether it runs successfully
        :return: None
        """
        if any(err in self.stderr.lower() for err in ("error", "fail")):
            raise ToolExecutionError('Command failed: {!r}\n stderr: {}'.format(self._command.command, self.stderr))
