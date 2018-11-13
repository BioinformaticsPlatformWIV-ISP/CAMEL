import os

from camel.app.error.invalidparametererror import InvalidParameterError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.samtools.samtools import Samtools


class SamtoolsSort(Samtools):
    """
    Sorts alignment files.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super(SamtoolsSort, self).__init__('samtools sort', '1.9', camel)

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        if 'BAM' not in self._tool_inputs:
            raise ValueError("No BAM input file found")
        super(Samtools, self)._check_input()

    def _check_parameters(self):
        """
        Checks the tool parameters.
        :return: None
        """
        if self._parameters['output_format'].value.upper() not in ('SAM', 'BAM'):
            raise InvalidParameterError("Invalid output format (BAM/SAM supported)")
        super(SamtoolsSort, self)._check_parameters()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()
        self._check_stderr()

    def __build_command(self):
        """
        Builds the command
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join(self._build_options()),
            self._tool_inputs['BAM'][0].path])

    def __set_output(self):
        """
        Sets the tool output.
        :return: None
        """
        output_path = os.path.join(self._folder, self._parameters['output_filename'].value)
        if not os.path.isfile(output_path):
            raise ToolExecutionError("Expected {} output not generated".format(self._name))
        output_key = self._parameters['output_format'].value.upper()
        self._tool_outputs[output_key] = [ToolIOFile(output_path)]

    def _check_command_output(self):
        """
        Checks if the command was executed successfully. Supersedes that of Tool class as samtools prints warnings to stderr.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Command execution failed (Exit code: {})".format(self._command.returncode))
