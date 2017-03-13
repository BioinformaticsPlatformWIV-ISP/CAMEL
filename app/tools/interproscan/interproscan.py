import shutil

import os

from app.tools.tool import Tool
from app.error.toolexecutionerror import ToolExecutionError
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.io.tooliofile import ToolIOFile
import logging


class Interproscan(Tool):
    """
    InterPro is a database which integrates together predictive information about proteins' function from a number of
    partner resources, giving an overview of the families that a protein belongs to and the domains and sites it
    contains. InterProScan is the software that can be used to query this database
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super(Interproscan, self).__init__('interproscan', '5.20-59.0', camel)
        self.__input_key = None

    def _execute_tool(self):
        """
        Runs InterProScan
        :return: None
        """
        self.__set_input_key()
        self.__build_command()
        self._execute_command()
        self.__set_output()
        shutil.rmtree(self.__get_temp_dir())

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA_Nucl or FASTA_Prot is required
        - Only one input file allowed per input key
        :return: None
        """
        super(Interproscan, self)._check_input()
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidInputSpecificationError('Invalid number of input keys given for InterProScan, only FASTA_Nucl '
                                                 'or FASTA_Prot allowed: {!r}'.format(self._tool_inputs))
        if 'FASTA_Nucl' not in self._tool_inputs and 'FASTA_Prot' not in self._tool_inputs:
            raise InvalidInputSpecificationError('No valid input key given for InterProScan, FASTA_Nucl or FASTA_Prot '
                                                 'needed: {!r}'.format(self._tool_inputs))
        for value in self._tool_inputs.values():
            if len(value) > 1:
                raise InvalidInputSpecificationError('Too many input files per key given for InterProScan '
                                                     '(max = 1): {!r}'.format(self._tool_inputs))

    def __set_input_key(self):
        """
        Sets the input key that is present in the inputs
        :return: None
        """
        self.__input_key = 'FASTA_Nucl' if 'FASTA_Nucl' in self._tool_inputs else 'FASTA_Prot'

    def __set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        if 'formats' in self._parameters:
            output_keys = self._parameters['formats'][1].split(',')
        else:
            output_keys = ['GFF3', 'XML', 'TSV']
        for key in output_keys:
            self._tool_outputs[key] = [ToolIOFile('{}.{}'.format(os.path.join(self._folder, self.__get_basename()), key.lower()))]

    def __get_basename(self):
        """
        Creates the basename for the output files
        :return: Basename
        """
        return self._tool_inputs[self.__input_key][0].basename

    def __build_input_string(self):
        """
        Creates the string with the input files
        :return: String with the input parameters
        """
        seqtype = 'n' if self.__input_key == 'FASTA_Nucl' else 'p'
        return '--input {} --seqtype {}'.format(self._tool_inputs[self.__input_key][0], seqtype)

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        input_string = self.__build_input_string()
        options_string = self.__build_options_string()
        self._command.command = ' '.join([self._tool_command, input_string, options_string])

    def __get_temp_dir(self):
        """
        Returns the path to the temporary directory and creates it if it does not yet exist
        :return: Path to temporary directory
        """
        if os.path.isfile(os.path.join(self._folder, 'temp')):
            raise IOError('A file with the name temp already exists!')
        if not os.path.isdir(os.path.join(self._folder, 'temp')):
            os.mkdir(os.path.join(self._folder, 'temp'))
        return os.path.join(self._folder, 'temp')

    def _check_command_output(self):
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if 'Unfortunately the web service has failed.' in self.stdout:
            logging.warning('The local lookup webservice was not reachable by InterProScan!')
        if self._command.returncode != 0:
            raise ToolExecutionError("Command execution failed for InterProScan (Exit code: {})".format(self._command.returncode))

    def __build_options_string(self):
        """
        Creates the string with all the specified parameters
        :return: String with command parameters
        """
        option_list = super(Interproscan, self)._build_options()
        option_list.append('--output-dir {}'.format(self._folder))
        option_list.append('--tempdir {}'.format(self.__get_temp_dir()))
        return '' if len(option_list) == 0 else ' '.join(option_list)
