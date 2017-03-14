import os

from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.error.toolexecutionerror import ToolExecutionError
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class Itsx(Tool):
    """
    ITSx is an open source software utility to extract the highly variable ITS1 and ITS2 subregions from ITS sequences,
    which is commonly used as a molecular barcode for e.g. fungi. As the inclusion of parts of the neighbouring, very
    conserved, ribosomal genes (SSU, 5S and LSU rRNA sequences) in the sequence identification process can lead to
    severely misleading results, ITSx identifies and extracts only the ITS regions themselves.
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super(Itsx, self).__init__('itsx', '1.0.11', camel)

    def _execute_tool(self):
        """
        Runs ITSx
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA key is required
        - Only one input file allowed
        - No other input keys are allowed
        :return: None
        """
        super(Itsx, self)._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('Not enough valid input files given for ITSx (FASTA required): {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidInputSpecificationError('Invalid number (max = 1) of FASTA files given for ITSx: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) > 2:
            raise InvalidInputSpecificationError('Too many input keys given voor ITSx (only FASTA allowed): {!r}'.format(self._tool_inputs))

    def __get_basename(self):
        """
        Returns the prefix that will be used in the output.
        :return: String with the prefix used in the output
        """
        infile = self._tool_inputs['FASTA'][0].basename
        return '{}_ITSx'.format(os.path.join(self._folder, os.path.splitext(infile)[0]))

    def __set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        basename = self.__get_basename()
        self._tool_outputs['TEXT_Summary'] = [ToolIOFile(basename + '.summary.txt')]
        self._tool_outputs['TEXT_NoDetection'] = [ToolIOFile(basename + '_no_detections.fasta')]
        self._tool_outputs['TEXT_Problematic'] = [ToolIOFile(basename + '.problematic.txt')]
        self._tool_outputs['TSV_Positions'] = [ToolIOFile(basename + '.positions.txt')]
        self._tool_outputs['GRAPH'] = [ToolIOFile(basename + '.graph')]
        self._tool_outputs['FASTA_Full'] = [ToolIOFile(basename + '.full.fasta')]
        self._tool_outputs['FASTA_ITS1'] = [ToolIOFile(basename + '.ITS1.fasta')]
        self._tool_outputs['FASTA_ITS2'] = [ToolIOFile(basename + '.ITS2.fasta')]
        self._tool_outputs['FASTA_NoDetection'] = [ToolIOFile(basename + '_no_detections.fasta')]
        if os.path.isfile(basename + '.chimeric.fasta'):
            self._tool_outputs['FASTA_Chimeric'] = [ToolIOFile(basename + '.chimeric.fasta')]

    def build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        return ' '.join(['-i {}'.format(self._tool_inputs['FASTA'][0]),
                         '-o {}'.format(self.__get_basename())])

    def __build_command(self):
        """
        Concatenates required parameters and options to build the command to run 'make.contigs'
        :return: None
        """
        input_string = self.build_input_string()
        options_string = ' '.join(self._build_options())
        self._command.command = '{} {} {}'.format(self._tool_command, input_string, options_string)

    def _check_command_output(self):
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Command execution for ITSx failed (Exit code: {})".format(self._command.returncode))
