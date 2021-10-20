from pathlib import Path

from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.samtools.samtools import Samtools


class SamtoolsIndex(Samtools):
    """
    Indexes sorted BAM files.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super().__init__('samtools index', '1.9', camel)

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        if 'BAM' not in self._tool_inputs:
            raise ValueError("No BAM input file found")
        if len(self._tool_inputs['BAM']) != 1:
            raise ValueError("Only one BAM input file is supported")
        super(Samtools, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        input_file_path = self.__symlink_input()
        self.__build_command(input_file_path)
        self._execute_command()
        self._check_stderr()
        self._tool_outputs['BAM'] = [ToolIOFile(Path(input_file_path))]

    def __symlink_input(self):
        """
        Create a symlink for the input. This avoids cluttering the directory of the input file. This can also avoid
        errors when there are no writing permissions on the directory of the input file.
        :return: Path to symlink input
        """
        if 'output_filename' in self._parameters:
            basename = self._parameters['output_filename'].value
        else:
            basename = self._tool_inputs['BAM'][0].basename
        new_path = self._folder / basename
        if (not new_path.is_symlink()) and (new_path != self._tool_inputs['BAM'][0].path):
            new_path.symlink_to(self._tool_inputs['BAM'][0].path)
        return new_path

    def __build_command(self, input_file_path):
        """
        Builds the command for this tool.
        :param input_file_path: Path to the input file
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join(self._build_options(excluded_parameters=['output_filename'])),
            str(input_file_path)])

    def _check_stderr(self):
        """
        Validates the stderr.
        :return: None
        """
        if 'unsorted positions' in self.stderr:
            raise ToolExecutionError('BAM file is not sorted.')
        super(SamtoolsIndex, self)._check_stderr()
