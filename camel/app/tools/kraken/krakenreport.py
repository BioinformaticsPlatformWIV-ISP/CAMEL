import os

from camel.app.tools.tool import Tool
from camel.app.io.tooliofile import ToolIOFile
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError


class KrakenReport(Tool):
    """
    Kraken is a system for assigning taxonomic labels to short DNA sequences, usually obtained through metagenomic
    studies. This class will generate a report of Kraken results across an entire sample.
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('kraken_report', '0.10.5', camel)

    def _execute_tool(self):
        """
        Runs Prinseq
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - TSV and DB keys are required
        - Only one input file allowed per key
        - No other input keys are allowed
        :return: None
        """
        if 'TSV' not in self._tool_inputs or 'DB' not in self._tool_inputs:
            raise InvalidInputSpecificationError('TSV or DB input keys are missing for Kraken-report: {!r}'.format(self._tool_inputs))
        for value in self._tool_inputs.values():
            if len(value) > 1:
                raise InvalidInputSpecificationError('More than one file per key given '
                                                     'for Kraken-report: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) > 2:
            raise InvalidInputSpecificationError('Too many input keys given for Kraken-report '
                                                 '(only TSV and DB allowed): {!r}'.format(self._tool_inputs))

    def __get_basename(self):
        """
        Returns the prefix that will be used in the output.
        :return: String with the prefix used in the output
        """
        infile = os.path.basename(self._tool_inputs['TSV'][0].path)
        return os.path.join(self._folder, os.path.splitext(infile)[0])

    def __set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs['TSV'] = [ToolIOFile(self.__get_basename() + '.report.tsv')]

    def __build_input_string(self):
        """
        Creates the string with the input and output files
        :return: String with the input parameters
        """
        command_parts = ['--db {}'.format(self._tool_inputs['DB'][0]),
                         '{}'.format(self._tool_inputs['TSV'][0]),
                         '> {}'.format(self.__get_basename() + '.report.tsv')]
        return ' '.join(command_parts)

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        options_string = ' '.join(self._build_options())
        self._command.command = '{} {} {}'.format(self._tool_command, self.__build_input_string(), options_string)
