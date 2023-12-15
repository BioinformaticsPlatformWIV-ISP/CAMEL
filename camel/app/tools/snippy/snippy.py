from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Snippy(Tool):

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('snippy', '4.6.0', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA_REF' not in self._tool_inputs:
            raise InvalidInputSpecificationError('Reference genome input is required.')
        if 'FASTA_CTGS' not in self._tool_inputs:
            raise InvalidInputSpecificationError('Fasta input is required')
        super()._check_input()

    def __build_command(self) -> None:
        """
        Builds the command line string.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            f"--reference {self._tool_inputs['FASTA_REF'][0].path}",
            f"--ctgs {self._tool_inputs['FASTA_CTGS'][0].path}",
            *self._build_options()]
        )

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__set_output()
        self.__build_command()
        self._execute_command()

    def __set_output(self) -> None:
        """
        Collects the output files of interest.
        :return: None
        """
        self._tool_outputs['BAM'] = [ToolIOFile(self.folder / self._parameters['output_directory'].value / 'snps.bam')]

    def _check_command_output(self) -> None:
        """
        Checks if the command was executed successfully.
        :return: None
        """
        # snippy logs to the stderr instead of the stdout
        if self._command.returncode != 0:
            raise ToolExecutionError("Command execution failed (Exit code: {})".format(self._command.returncode))
