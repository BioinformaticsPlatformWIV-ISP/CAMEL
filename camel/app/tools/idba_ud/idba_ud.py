import os
from Bio import SeqIO

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class IdbaUd(Tool):
    """
    IDBA is the basic iterative de Bruijn graph assembler for second-generation sequencing reads. IDBA-UD, an extension
    of IDBA, is designed to utilize paired-end reads to assemble low-depth regions and use progressive depth on contigs
    to reduce errors in high-depth regions.
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('idba_ud', '1.1.1', camel)

    def _execute_tool(self):
        """
        Runs IDBA_UD
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - Only FASTA key allowed
        - Up to 5 input files allowed (up to 5th level scaffolds)
        :return: None
        """
        super(IdbaUd, self)._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('Invalid input key given for IDBA_UD, FASTA is required: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidInputSpecificationError('Invalid number of input keys given for IDBA_UD, only FASTA is allowed: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs['FASTA']) > 5:
            raise InvalidInputSpecificationError('Invalid number of files given for IDBA_UD, maximum is 5: {!r}'.format(self._tool_inputs))

    def __set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs['FASTA_Contig'] = [ToolIOFile(os.path.join(self._folder, 'contig.fa'))]
        self._tool_outputs['FASTA_Scaffold'] = [ToolIOFile(os.path.join(self._folder, 'scaffold.fa'))]
        self._tool_outputs['LOG'] = [ToolIOFile(os.path.join(self._folder, 'log'))]

    def __build_input_string(self):
        """
        Creates the string with the input files
        :return: String with the input parameters
        """
        parts = []
        if len(max(SeqIO.parse(self._tool_inputs['FASTA'][0].path, 'fasta'), key=len)) > 302:
            parts.append('-l {}'.format(self._tool_inputs['FASTA'][0]))
        else:
            parts.append('-r {}'.format(self._tool_inputs['FASTA'][0]))
        for i in range(1, len(self._tool_inputs['FASTA'])):
            parts.append('--read_level_{} {}'.format(i+1, self._tool_inputs['FASTA'][i]))
        return ' '.join(parts)

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        input_string = self.__build_input_string()
        out_string = '-o {}'.format(os.path.join(self._folder))
        options_string = ' '.join(self._build_options())
        self._command.command = ' '.join([self._tool_command, input_string, out_string, options_string])

    def _check_command_output(self):
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Command execution failed for IDBA_UD (Exit code: {})".format(self._command.returncode))
