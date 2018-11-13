import os

import re

from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.samtools.samtools import Samtools


class SamtoolsFlagstat(Samtools):
    """
    Calculates Simple BAM/SAM file statistics.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super(SamtoolsFlagstat, self).__init__('samtools flagstat', '1.9', camel)

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        if 'BAM' not in self._tool_inputs:
            raise ValueError("No BAM input file found")
        super(Samtools, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_informs()
        self.__set_output()
        self._check_stderr()

    def __build_command(self):
        """
        Builds the command
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            self._tool_inputs['BAM'][0].path])

    def __set_informs(self):
        """
        Sets the informs for this tool.
        :return: None
        """
        lines = self.stdout.split('\n')
        self._informs['total'] = SamtoolsFlagstat.__parse_output_line(lines[0])
        self._informs['secondary'] = SamtoolsFlagstat.__parse_output_line(lines[1])
        self._informs['supplementary'] = SamtoolsFlagstat.__parse_output_line(lines[2])
        self._informs['duplicates'] = SamtoolsFlagstat.__parse_output_line(lines[3])
        self._informs['mapped'] = SamtoolsFlagstat.__parse_output_line(lines[4])
        self._informs['paired'] = SamtoolsFlagstat.__parse_output_line(lines[5])
        self._informs['read1'] = SamtoolsFlagstat.__parse_output_line(lines[6])
        self._informs['read2'] = SamtoolsFlagstat.__parse_output_line(lines[7])
        self._informs['properly_paired'] = SamtoolsFlagstat.__parse_output_line(lines[8])
        self._informs['singletons'] = SamtoolsFlagstat.__parse_output_line(lines[10])

    @staticmethod
    def __parse_output_line(line):
        """
        Parses a line of flagstat output
        :param line: Flagstat output line
        :return: Line values
        """
        m = re.match('^(\d+) \+ (\d+).*', line)
        if m is None:
            raise ValueError("Cannot parse: '{}'".format(line))
        return int(m.group(1)), int(m.group(2))

    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """
        if 'output_filename' in self._parameters:
            output_path = os.path.join(self._folder, self._parameters['output_filename'].value)
            try:
                with open(output_path, 'w') as handle:
                    handle.write(self.stdout)
                self._tool_outputs['TXT'] = [ToolIOFile(output_path)]
            except IOError:
                raise ToolExecutionError("Cannot create output file")
