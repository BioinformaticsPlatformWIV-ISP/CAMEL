import os

from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class Fqextract(Tool):
    """
    Extract reads from FASTQ files.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super(Fqextract, self).__init__('fqextract', '1.1', camel)

    def _check_input(self):
        """
        Checks the tool input.
        :return: None
        """
        if 'FASTQ' not in self._tool_inputs:
            raise ValueError("No FASTQ input found.")
        if 'TXT' not in self._tool_inputs:
            raise ValueError("No sequence input list found")
        if len(self._tool_inputs['FASTQ']) != 1 or len(self._tool_inputs['TXT']) != 1:
            raise ValueError("Only one input is supported.")
        super(Fqextract, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __build_command(self):
        """
        Builds the command
        :return: None
        """
        self._command.command = ' '.join(
            [self._tool_command,
             '-l {}'.format(self._tool_inputs['TXT'][0].path),
             ' '.join(self._build_options()),
             self._tool_inputs['FASTQ'][0].path])

    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['FASTQ'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output'].value))]
