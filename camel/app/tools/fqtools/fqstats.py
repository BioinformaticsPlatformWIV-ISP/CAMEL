import re

from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.tools.tool import Tool


class Fqstats(Tool):
    """
    Reports the number of sequences & the number of bases from FASTQ files.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super(Fqstats, self).__init__('fqstats', '1.1',  camel)

    def _check_input(self):
        """
        Checks the tool input.
        :return: None
        """
        if 'FASTQ' not in self._tool_inputs:
            raise ValueError("No FASTQ input found.")
        super(Fqstats, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        for fastq_file in self._tool_inputs['FASTQ']:
            self.__build_command(fastq_file)
            self._execute_command()
            self.__add_informs(fastq_file)

    def __build_command(self, input_file):
        """
        Builds the command
        :param input_file: Input FASTQ file
        :return: None
        """
        self._command.command = '{} {}'.format(self._tool_command, input_file.path)

    def __add_informs(self, input_file):
        """
        Adds the informs.
        :param input_file: Input file
        :return: None
        """
        m = re.match('Found (\d+) sequences, (\d+) bases', self.stdout.strip())
        if not m:
            raise ToolExecutionError("Error retrieving fqstats from {}.".format(input_file.basename))
        self._informs[input_file.path] = {'nb_of_sequences': int(m.group(1)), 'nb_of_bases': int(m.group(2))}
