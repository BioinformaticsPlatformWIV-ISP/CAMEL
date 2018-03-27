import os

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class BcftoolsView(Tool):
    """
    VCF/BCF conversion, view, subset and filter VCF/BCF files.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(BcftoolsView, self).__init__('bcftools view', '1.3.1', camel)

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        if not any(key in self._tool_inputs for key in ('BCF', 'VCF', 'VCF_GZ')):
            raise InvalidInputSpecificationError("No input file found (BCF / VCF_GZ / VCF supported)")
        super(BcftoolsView, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __get_input_file_path(self):
        """
        Returns the path to the input file.
        :return: Input file path
        """
        if 'VCF_GZ' in self._tool_inputs:
            return self._tool_inputs['VCF_GZ'][0].path
        elif 'VCF' in self._tool_inputs:
            return self._tool_inputs['VCF'][0].path
        else:
            return self._tool_inputs['BCF'][0].path

    def __get_output_format_option(self):
        """
        Returns the output command line option.
        :return: Command line option
        """
        if self._parameters['output_format'].value == 'VCF':
            return '--output-type z' if 'compress_output' in self._parameters else '--output-type v'
        else:
            return '--output-type b' if 'compress_output' in self._parameters else '--output-type u'

    def __get_output_key(self):
        """
        Returns the output key.
        :return: Output key
        """
        if self._parameters['output_format'].value == 'VCF':
            return 'VCF_GZ' if 'compress_output' in self._parameters else 'VCF'
        else:
            return 'BCF'

    def __build_command(self):
        """
        Builds the command.
        :return: None
        """
        command_parts = [
            self._tool_command,
            self.__get_input_file_path(),
            self.__get_output_format_option()
        ]
        command_parts += self._build_options(['output_format', 'compress_output'])
        self._command.command = ' '.join(command_parts)

    def _check_command_output(self):
        """
        Checks if the command executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Error executing bcftools view: {}".format(self.stderr))

    def __set_output(self):
        """
        Sets the tool output.
        :return: None
        """
        self._tool_outputs[self.__get_output_key()] = [
            ToolIOFile(os.path.join(self._folder, self._parameters['output_filename'].value))]
