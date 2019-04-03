import os

from camel.app.error.invalidparametererror import InvalidParameterError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.samtools.samtools import Samtools


class SamtoolsView(Samtools):
    """
    SAM <-> BAM Conversion
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super().__init__('samtools view', '1.9', camel)
        self.__input_key = None

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        if 'SAM' in self._tool_inputs:
            self.__input_key = 'SAM'
        elif 'BAM' in self._tool_inputs:
            self.__input_key = 'BAM'
        else:
            raise ValueError("No input file found")
        super(Samtools, self)._check_input()

    def _check_parameters(self):
        """
        Checks the tool parameters.
        :return: None
        """
        if self._parameters['output_format'].value.upper() not in ('SAM', 'BAM'):
            raise InvalidParameterError("Invalid output format (BAM/SAM supported)")
        super(SamtoolsView, self)._check_parameters()

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
        Builds the command for this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join(self._build_options(['regions'])),
            self._tool_inputs[self.__input_key][0].path])
        if 'regions' in self._parameters:
            self._command.command += ' "{}"'.format(self._parameters['regions'].value)

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

    def _check_stderr(self):
        """
        Validates the stderr.
        :return: None
        """
        if 'only works for indexed' in self.stderr:
            raise ToolExecutionError("Can only extract regions from indexed BAM files")
        super(SamtoolsView, self)._check_stderr()
