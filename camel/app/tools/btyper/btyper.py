from camel.app.camel import Camel
from camel.app.tools.tool import Tool

class Btyper(Tool):

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Btyper', '3.2.0', camel)

    def _check_input(self) -> None:
        """
        Checks whether the required input files are specified.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise ValueError('No FASTA input found')
        super(Btyper, self)._check_input()

    def _build_command(self) -> None:
        """
        Build the command to run tool
        :return: None
        """
        self._command.command = f'{self._tool_command} ' \
                                f'--input {self._tool_inputs["FASTA"][0]} ' \
                                f'{" ".join(self._build_options())}'

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._check_input()
        self._build_command()
        self._execute_command()
