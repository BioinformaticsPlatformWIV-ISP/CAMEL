import os

import re

from camel.app.error.invalidparametererror import InvalidParameterError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Fqconvert(Tool):
    """
    Converts fastq files scoring system.
    """

    def __init__(self, camel):
        """
        Initialize FastQC tool
        :param camel: Camel instance
        """
        self._mode = None
        super(Fqconvert, self).__init__('fqconvert', '1.1', camel)

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self._add_informs()
        self._add_output_files()

    def _check_parameters(self):
        """
        Checks the parameters.
        :return: None
        """
        if 'detect_only' in self._parameters:
            self._mode = 'detection'
        else:
            self._mode = 'convert'
            if 'output' not in self._parameters:
                raise InvalidParameterError("Output parameter is required when doing conversion")
        super(Fqconvert, self)._check_parameters()

    def _check_input(self):
        """
        Checks the tool input.
        :return: None
        """
        if 'FASTQ' not in self._tool_inputs:
            raise ValueError("No FASTQ input found")
        if self._mode == 'convert' and len(self._tool_inputs['FASTQ']) > 1:
            raise ValueError("Only one input is supported when doing conversion")
        super(Fqconvert, self)._check_input()

    def __build_command(self):
        """
        Builds the command.
        :return: Command
        """
        self._command.command = ' '.join(
            [self._tool_command,
             ' '.join(self._build_options()),
             ' '.join([i.path for i in self._tool_inputs['FASTQ']])])

    def _check_command_output(self):
        """
        Checks the command output.
        :return: None
        """
        for line in self.stderr.splitlines():
            if not re.match('.*: Found.*', line):
                raise ToolExecutionError("Error executing fqconvert: {}".format(line.strip()))

    def _add_informs(self):
        """
        Adds the informs.
        :return: None
        """
        for line in self.stderr.splitlines():
            m = re.match('^(.+): Found (.+)$', line.strip())
            self._informs[m.group(1)] = m.group(2)

    def _add_output_files(self):
        """
        Adds the output files.
        :return: None
        """
        if self._mode == 'convert':
            self._tool_outputs['FASTQ'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output'].value))]
