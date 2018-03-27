import os

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Deconseq(Tool):
    """
    The DeconSeq tool can be used to automatically detect and efficiently remove sequence contaminations from genomic
    and metagenomic datasets. It is easily configurable and provides a user-friendly interface.
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super(Deconseq, self).__init__('deconseq', '0.4.3', camel)
        self._input_key = None

    def _execute_tool(self):
        """
        Runs Deconseq
        :return: None
        """
        self.__set_input_key()
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTQ or FASTA key is required
        - Only one input file allowed
        - No other input keys are allowed
        :return: None
        """
        for key, value in self._tool_inputs.items():
            if key not in ['FASTQ', 'FASTA']:
                raise InvalidInputSpecificationError('Illegal input key given for DeconSeq, '
                                                     'only FASTQ or FASTA allowed: {!r}'.format(self._tool_inputs))
            if len(value) != 1:
                raise InvalidInputSpecificationError('Illegal number of input files (max = 1) '
                                                     'provided for DeconSeq: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidInputSpecificationError('Too many input keys given for DeconSeq, '
                                                 'only FASTQ or FASTA allowed: {!r}'.format(self._tool_inputs))

    def __get_basename(self):
        """
        Returns the prefix that will be used in the output.
        :return: String with the prefix used in the output
        """
        infile = os.path.basename(self._tool_inputs[self._input_key][0].path)
        return os.path.join(self._folder, os.path.splitext(infile)[0]) + '.deconseq'

    def __get_extension(self):
        """
        Returns the extension that the output file will have
        :return: Extension of the output file
        """
        if self._input_key == 'FASTQ':
            return '.fq'
        else:
            return '.fa'

    def __set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        basename = self.__get_basename()
        extension = self.__get_extension()
        self._tool_outputs[self._input_key + '_Clean'] = [ToolIOFile(basename + '_clean' + extension)]
        self._tool_outputs[self._input_key + '_Cont'] = [ToolIOFile(basename + '_cont' + extension)]
        if 'dbs_retain' in self._parameters:
            self._tool_outputs[self._input_key + '_Both'] = [ToolIOFile(basename + '_both' + extension)]

    def __build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        command_parts = ['-f {}'.format(self._tool_inputs[self._input_key][0]),
                         '-id {}'.format(os.path.basename(self.__get_basename())),
                         '-out_dir {}'.format(self._folder)]
        return ' '.join(command_parts)

    def __set_input_key(self):
        """
        Sets the instance variable self._input_key
        :return: None
        """
        self._input_key = list(self._tool_inputs.keys())[0]

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        input_string = self.__build_input_string()
        options_string = ' '.join(self._build_options())
        self._command.command = '{} {} {}'.format(self._tool_command, input_string, options_string)

    def _check_command_output(self):
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Command execution failed (Exit code: {})".format(self._command.returncode))
