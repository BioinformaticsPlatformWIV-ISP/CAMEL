import logging
import os

from app.error.invalidparametererror import InvalidParameterError
from app.error.toolexecutionerror import ToolExecutionError
from app.io.tooliofile import ToolIOFile
from app.tools.samtools.samtools import Samtools


class SamtoolsFastaIndex(Samtools):
    """
    Indexes FASTA files.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super(SamtoolsFastaIndex, self).__init__('samtools faidx', '1.3', camel)

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise ValueError("No FASTA input file found")
        if len(self._tool_inputs['FASTA']) != 1:
            raise ValueError("Only one FASTA input file is supported.")
        super(Samtools, self)._check_input()

    def _check_parameters(self):
        """
        Checks the parameters.
        :return: None
        """
        if 'regions' in self._parameters and 'output_filename' not in self._parameters:
            raise InvalidParameterError("Cannot extract regions without output filename")
        super(SamtoolsFastaIndex, self)._check_parameters()

    def __symlink_input(self):
        """
        Creates a symlink for the input.
        :return: Path to the symlink of the input
        """
        symlink_location = os.path.join(self._folder, self._tool_inputs['FASTA'][0].basename)
        try:
            os.symlink(self._tool_inputs['FASTA'][0].path, symlink_location)
        except OSError:
            pass
        return symlink_location

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        fasta_file = self.__symlink_input()
        self.__build_command(fasta_file)
        self._execute_command()
        self._check_stderr()
        if 'regions' in self._parameters:
            self._tool_outputs['FASTA'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output_filename'].value))]
        else:
            self._tool_outputs['FASTA'] = [ToolIOFile(fasta_file)]

    def __build_command(self, fasta_file):
        """
        Builds the command for this tool.
        :param fasta_file: FASTA file
        :return: None
        """
        self._command.command = ' '.join([self._tool_command, fasta_file])
        if 'output_filename' in self._parameters and 'regions' in self._parameters:
            logging.info("Extracting regions from FASTA file, file should already be indexed.")
            self._command.command += ' {} > {}'.format(
                self._parameters['regions'].value, self._parameters['output_filename'].value)

    def _check_stderr(self):
        """
        Checks the command stderr output.
        :return: None
        """
        if 'build FASTA index' in self.stderr:
            raise ToolExecutionError("Cannot extract regions from an unindexed FASTA file.")
