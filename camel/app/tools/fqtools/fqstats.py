import re
from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.tools.tool import Tool


class FqStats(Tool):
    """
    Reports the number of sequences & the number of bases from input FASTQ files.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super().__init__('fqstats', '1.1',  camel)

    def _check_input(self) -> None:
        """
        Check if the provided tool input is valid.
        :return: None
        """
        if 'FASTQ' not in self._tool_inputs:
            raise InvalidInputSpecificationError("FASTQ input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._informs['stats'] = []
        for fastq_file in self._tool_inputs['FASTQ']:
            self.__build_command(fastq_file.path)
            self._execute_command()
            self.__add_informs(fastq_file. path)

    def __build_command(self, fq_in: Path) -> None:
        """
        Builds the command
        :param fq_in: Input FASTQ file
        :return: None
        """
        self._command.command = ' '.join([self._tool_command, str(fq_in)])

    def __add_informs(self, fq_in: Path) -> None:
        """
        Adds the informs.
        :param fq_in: Input FASTQ file
        :return: None
        """
        m = re.match(r'Found (\d+) sequences, (\d+) bases', self.stdout.strip())
        if not m:
            raise ToolExecutionError(f"Error parsing fqstats for file: {fq_in.name}.")
        self._informs['stats'].append({
            'path': str(fq_in),
            'nb_of_sequences': int(m.group(1)),
            'nb_of_bases': int(m.group(2))
        })
