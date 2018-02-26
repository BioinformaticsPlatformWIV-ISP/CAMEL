import os

import logging

from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class BlastFormatterLoop(Tool):
    """
    Formats BLAST output.
    """

    def __init__(self, camel):
        """
        Initialize tool.
        :param camel: Camel instance
        :return: None
        """
        super(BlastFormatterLoop, self).__init__('blast_formatter (looping)', '2.6.0', camel)

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        if len(self._tool_inputs) == 0:
            logging.info("Not enough inputs for {}".format(self.name))
            return
        output_key = self.__get_output_key()
        self._tool_outputs[output_key] = []

        for blast_archive in self._tool_inputs['ASN']:
            output_name = self.__get_output_name(blast_archive)
            self.__run_blast_formatter(blast_archive, output_name)
            self._tool_outputs[output_key].append(ToolIOFile(os.path.join(self._folder, output_name)))

    def _check_input(self):
        """
        Check if the required input files were specified.
        :return: None
        """
        if len(self._tool_inputs) == 0:
            pass
        elif 'ASN' not in self._tool_inputs:
            raise ValueError("No BLAST archive input found")
        super(BlastFormatterLoop, self)._check_input()

    def __check_command_output(self):
        """
        Checks the command stdout and stderr output for errors.
        :return: None
        """
        if 'error' in self._command.stderr.lower():
            raise Exception("Problem running blast_formatter: {}".format(self._command.stderr))

    def __run_blast_formatter(self, blast_archive, output_name):
        """
        Runs this tool once with the given parameters.
        :param blast_archive: BLAST archive to format
        :param output_name: Name of the output file
        :return: None
        """
        self._command.command = '{} -archive {} -out {} {}'.format(
            self._tool_command, blast_archive.path, output_name, ' '.join(self._build_options()))
        self._execute_command()
        self.__check_command_output()

    def __get_output_name(self, input_file):
        """
        Generates the default output name.
        :return: Output name
        """
        base_filename = os.path.splitext(input_file.basename)[0]
        return '{}_formatted.{}'.format(base_filename, self.__get_output_key().lower())

    def __get_output_key(self):
        """
        Returns the output key.
        :return:
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

    def _check_command_output(self):
        """
        Checks the command output for errors.
        :return: None
        """
        if 'error' in self.stderr.lower():
            raise ValueError("Error executing {}: {}".format(self.name, self._command.stderr.strip()))
