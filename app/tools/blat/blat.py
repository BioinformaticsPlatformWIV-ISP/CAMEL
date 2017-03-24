import os

from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.error.toolexecutionerror import ToolExecutionError
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class Blat(Tool):
    """
    BLAT is a bioinformatics software a tool which performs rapid mRNA/DNA and cross-species protein alignments. BLAT
    is more accurate and 500 times faster than popular existing tools for mRNA/DNA alignments and 50 times faster for
    protein alignments at sensitivity settings typically used when comparing vertebrate sequences.
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super(Blat, self).__init__('blat', '36x2', camel)

    def _execute_tool(self):
        """
        Runs BLAT
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - DB_DNA, DB_Prot or DB_DNAX required
        - FASTA_DNA, FASTA_Prot, FASTA_RNA, FASTA_DNAX, or FASTA_RNAX required
        - OOC optional
        - Only one input file allowed per input key
        :return: None
        """
        super(Blat, self)._check_input()
        if not 2 <= len(self._tool_inputs.keys()) <= 3:
            raise InvalidInputSpecificationError('Invalid number (min 2, max 3) of input keys given for BLAT: {!r}'.format(self._tool_inputs))
        if not any(key in self._tool_inputs for key in ['DB_DNA', 'DB_Prot', 'DB_DNAX']):
            raise InvalidInputSpecificationError('No database given for BLAT: {!r}'.format(self._tool_inputs))
        if not any(key in self._tool_inputs for key in ['FASTA_DNA', 'FASTA_Prot', 'FASTA_RNA', 'FASTA_DNAX', 'FASTA_RNAX']):
            raise InvalidInputSpecificationError('No query file given for BLAT: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) == 3 and 'OOC' not in self._tool_inputs:
            raise InvalidInputSpecificationError('Invalid input key given for BLAT: {!r}'.format(self._tool_inputs))
        for value in self._tool_inputs.values():
            if len(value) > 1:
                raise InvalidInputSpecificationError('Too many input files per key (max = 1) given for BLAT: {!r}'.format(self._tool_inputs))

    def __get_input_key(self):
        """
        Returns the input key that is present in the inputs
        :return: Input key
        """
        return [x for x in self._tool_inputs if x in ['FASTA_DNA', 'FASTA_Prot', 'FASTA_RNA', 'FASTA_DNAX', 'FASTA_RNAX']][0]

    def __get_db_key(self):
        """
        Returns the database key key that is present in the inputs
        :return: Database key
        """
        return [x for x in self._tool_inputs if x in ['DB_DNA', 'DB_Prot', 'DB_DNAX']][0]

    def __set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs[self.__get_output_key()] = [ToolIOFile(os.path.join(self._folder, 'output{}'.format(self.__get_output_extension())))]

    def __get_output_key(self):
        """
        Returns the output key that has to be used based on the specified parameters in the database
        :return: Output key to use
        """
        try:
            return 'TSV' if self._parameters['out'].value in ['psl', 'pslx', 'blast', 'blast8', 'blast9'] else \
                self._parameters['out'].value.upper()
        except KeyError:
            return 'TSV'

    def __get_output_extension(self):
        """
        Returns the extension of the output file
        :return: Output file extension
        """
        return '.{}'.format(self._parameters['out'].value) if 'out' in self._parameters else '.psl'

    def __build_input_string(self):
        """
        Creates the string with the input files
        :return: String with the input parameters
        """
        return '{} {}'.format(self._tool_inputs[self.__get_db_key()][0], self._tool_inputs[self.__get_input_key()][0])

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        input_string = self.__build_input_string()
        options_string = ' '.join(self._build_options(delimiter=''))
        out_string = 'output{}'.format(self.__get_output_extension())
        self._command.command = ' '.join([self._tool_command, input_string, options_string, out_string])

    def _check_command_output(self):
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Command execution failed (Exit code: {})".format(self._command.returncode))
