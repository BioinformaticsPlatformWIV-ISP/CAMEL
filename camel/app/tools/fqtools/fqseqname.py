from camel.app.tools.tool import Tool


class Fqseqname(Tool):
    """
    Extract sequences names for FASTQ files.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super().__init__('fqseqname', '1.1', camel)

    def _check_input(self):
        """
        Checks the tool input.
        :return: None
        """
        if 'FASTQ' not in self._tool_inputs:
            raise ValueError("No FASTQ input found.")
        super(Fqseqname, self)._check_input()

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
        self._informs[input_file.path] = [line.strip() for line in self.stdout.splitlines()]
