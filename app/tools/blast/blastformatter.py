import os

from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class BlastFormatter(Tool):
    """
    Formats BLAST output.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(BlastFormatter, self).__init__('blast_formatter', '2.6.0', camel)

    def _check_input(self):
        """
        Checks whether the required input files are specified.
        :return: None
        """
        if 'ASN' not in self._tool_inputs:
            raise ValueError('No blast archive input found')
        super(BlastFormatter, self)._check_input()

    def _execute_tool(self):
        """
        Runs Blast.
        :return: None
        """
        self.__output_key = self.__get_output_key()
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __get_output_key(self):
        """
        Returns the output key.
        :return: Output key
        """
        if 'create_html_output' in self._parameters:
            return 'HTML'
        output_format = self._parameters['output_format'].value
        if output_format == '5':
            return 'XML'
        elif '6' in output_format or '7' in output_format:
            return 'TSV'
        elif output_format in ('8', '9', '11'):
            return 'ASN'
        elif '10' in output_format:
            return 'CSV'
        elif output_format == '12':
            return 'JSON'
        else:
            return 'TXT'

    def __build_command(self):
        """
        Builds the command line string.
        :return: None
        """
        blast_archive = self._tool_inputs['ASN'][0].path
        output_name = self.__get_output_name()
        self._command.command = '{} -archive {} -out {} {}'.format(
            self._tool_command, blast_archive, output_name, ' '.join(self._build_options()))

    def __get_output_name(self):
        """
        Generates the default output name.
        :return: Output name
        """
        base_filename = os.path.splitext(self._tool_inputs['ASN'][0].basename)[0]
        return '{}_formatted.{}'.format(base_filename, self.__output_key.lower())

    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs[self.__output_key] = [ToolIOFile(os.path.join(self._folder, self.__get_output_name()))]

    def _check_command_output(self):
        """
        Checks the command output for errors.
        :return: None
        """
        if 'error' in self._command.stderr.lower():
            raise ValueError("Error executing {}: {}".format(self.name, self._command.stderr.strip()))
