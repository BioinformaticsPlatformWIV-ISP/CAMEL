from app.error.toolexecutionerror import ToolExecutionError
from app.tools.samtools.samtools import Samtools


class SamtoolsIndex(Samtools):
    """
    Indexes sorted BAM files.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super(SamtoolsIndex, self).__init__('samtools index', '1.3', camel)

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        if 'BAM' not in self._tool_inputs:
            raise ValueError("No BAM input file found")
        super(Samtools, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self._check_stderr()

    def __build_command(self):
        """
        Builds the command for this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join(self._build_options()),
            self._tool_inputs['BAM'][0].path])

    def _check_stderr(self):
        """
        Validates the stderr.
        :return: None
        """
        if 'unsorted positions' in self.stderr:
            raise ToolExecutionError('BAM file is not sorted.')
        super(SamtoolsIndex, self)._check_stderr()
