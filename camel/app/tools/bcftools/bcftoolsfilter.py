import os

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class BcftoolsFilter(Tool):
    """
    Filtering of VCF/BCF files using fixed thresholds.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('bcftools filter', '1.3.1', camel)

    def _check_input(self):
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'VCF_GZ' not in self._tool_inputs:
            raise InvalidInputSpecificationError("VCF_GZ input is required")
        super()._check_input()

    def __get_output_key(self):
        """
        Returns the output key.
        :return: Output key
        """
        if 'output_type' not in self._parameters:
            return 'VCF'
        output_type = self._parameters['output_type'].value
        if output_type == 'b':
            return 'BCF_GZ'
        elif output_type == 'u':
            return 'BCF'
        elif output_type == 'z':
            return 'VCF_GZ'
        else:
            return 'VCF'

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self._build_command()
        self._execute_command()
        self._tool_outputs[self.__get_output_key()] = [ToolIOFile(os.path.join(
            self._folder, self._parameters['output_filename'].value))]

    def _build_command(self):
        """
        Builds the command that is called.
        :return: None
        """
        command_parts = [self._tool_command]
        if 'BED' in self._tool_inputs:
            command_parts.append('--targets-file {}'.format(self._tool_inputs['BED'][0].path))
        command_parts.extend(self._build_options())
        command_parts.append(self._tool_inputs['VCF_GZ'][0].path)
        self._command.command = ' '.join(command_parts)
