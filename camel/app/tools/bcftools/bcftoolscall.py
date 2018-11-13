import logging
import os

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.invalidparametererror import InvalidParameterError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class BcftoolsCall(Tool):
    """
    SNP/indel variant calling from VCF/BCF. To be used in conjunction with samtools mpileup.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(BcftoolsCall, self).__init__('bcftools call', '1.9', camel)

    def _check_parameters(self):
        """
        Checks the parameters.
        :return: None
        """
        if self._parameters['calling_method'].value not in ('consensus', 'multiallelic'):
            raise InvalidParameterError("Unrecognized snp calling method: {}".format(
                self._parameters['calling_method'].value))
        if 'ploidy' not in self._parameters:
            logging.warning("Ploidy not specified will assume all sites are diploid.")

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        if not any(key in self._tool_inputs for key in ('BCF', 'VCF_GZ')):
            raise InvalidInputSpecificationError("No input file found (BCF / VCF_GZ supported)")
        super(BcftoolsCall, self)._check_input()

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
        else:
            return self._tool_inputs['BCF'][0].path

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
        if self._parameters['calling_method'].value == 'consensus':
            command_parts.append('--consensus-caller')
        else:
            command_parts.append('--multiallelic-caller')

        if 'TXT_RG' in self._tool_inputs:
            command_parts.append('--regions-file {}'.format(self._tool_inputs['TXT_RG'][0].path))
        if 'TXT_SAMPLES' in self._tool_inputs:
            command_parts.append('--samples-file {}'.format(self._tool_inputs['TXT_SAMPLES'][0].path))

        command_parts += self._build_options(['calling_method', 'output_format', 'compress_output'])
        self._command.command = ' '.join(command_parts)

    def _check_command_output(self):
        """
        Checks if the command executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Error executing bcftools call: {}".format(self.stderr))

    def __get_output_format_option(self):
        """
        Returns the output command line option.
        :return: Command line option
        """
        if self._parameters['output_format'].value == 'VCF':
            return '-O z' if 'compress_output' in self._parameters else '-O v'
        else:
            return '-O b' if 'compress_output' in self._parameters else '-O u'

    def __get_output_key(self):
        """
        Returns the output key.
        :return: Output key
        """
        if self._parameters['output_format'].value == 'VCF':
            return 'VCF_GZ' if 'compress_output' in self._parameters else 'VCF'
        else:
            return 'BCF'

    def __set_output(self):
        """
        Sets the tool output.
        :return: None
        """
        self._tool_outputs[self.__get_output_key()] = [
            ToolIOFile(os.path.join(self._folder, self._parameters['output_filename'].value))]
