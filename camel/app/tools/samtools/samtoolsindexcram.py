import os

from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.samtools.samtools import Samtools


class SamtoolsIndexCram(Samtools):
    """
    Indexes sorted CRAM files.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super().__init__('samtools index', '1.9', camel)

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        if 'CRAM' not in self._tool_inputs:
            raise ValueError("No CRAM input file found")
        if len(self._tool_inputs['CRAM']) != 1:
            raise ValueError("Only one CRAM input file is supported")

        if 'FASTA_REF' not in self._tool_inputs:
            raise ValueError("No FASTA_REF input file found")

        super(Samtools, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        input_file_path = self.__symlink_input()
        self.__build_command(input_file_path)
        self._execute_command()
        self._check_stderr()
        self._tool_outputs['CRAI'] = [ToolIOFile(input_file_path)]

    def __symlink_input(self):
        """
        Create a symlink for the input. This avoids cluttering the directory of the input file. This can also avoid
        errors when there are no writing permissions on the directory of the input file.
        :return: Path to symlink input
        """
        if 'output_filename' in self._parameters:
            basename = self._parameters['output_filename'].value
        else:
            basename = self._tool_inputs['CRAM'][0].basename
        new_path = os.path.join(self._folder, basename)
        if (not os.path.islink(new_path)) and (new_path != self._tool_inputs['CRAM'][0].path):
            os.symlink(self._tool_inputs['CRAM'][0].path, new_path)
        return new_path

    def __build_command(self, input_file_path):
        """
        Builds the command for this tool.
        seq_cache_populate.pl: Create REF_CACHE. Used when indexing a CRAM
        :param input_file_path: Path to the input file
        :return: None
        """
        self._command.command = ' '.join([
            "seq_cache_populate.pl -root ./ref/cache " + self._tool_inputs['FASTA_REF'][0].path + "; ",
            "export REF_PATH=:; ",
            "export REF_CACHE=./ref/cache/%2s/%2s/%s; ",
            self._tool_command,
            ' '.join(self._build_options(excluded_parameters=['output_filename'])),
            input_file_path])

    def _check_stderr(self):
        """
        Validates the stderr.
        :return: None
        """
        if 'unsorted positions' in self.stderr:
            raise ToolExecutionError('BAM file is not sorted.')
        super(SamtoolsIndexCram, self)._check_stderr()
