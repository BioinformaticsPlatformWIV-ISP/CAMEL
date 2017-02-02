import os

from app.error.toolexecutionerror import ToolExecutionError
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
        super(Samtools, self)._check_input()

    def __symlink_input(self):
        """
        Creates a symlink for the input.
        :return: Path to the symlink of the input
        """
        symlink_location = os.path.join(self._folder, self._tool_inputs['FASTA'][0].basename)
        os.symlink(self._tool_inputs['FASTA'][0].path, symlink_location)
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

    def __build_command(self, fasta_file):
        """
        Builds the command for this tool.
        :param fasta_file: FASTA file
        :return: None
        """
        self._command.command = ' '.join([self._tool_command, fasta_file])
        if 'output' in self._parameters and 'regions' in self._parameters:
            self._command.command += ' {} > {}'.format(
                self._parameters['regions'].value, self._parameters['output'].value)

    def _check_stderr(self):
        """
        Checks the command stderr output.
        :return: None
        """
        if 'build FASTA index' in self.stderr:
            raise ToolExecutionError("Cannot extract regions from an unindexed FASTA file.")
