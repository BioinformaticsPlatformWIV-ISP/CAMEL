from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class ALE2Wiggle(Tool):

    """
    ALE is the Assembly Likelihood Evaluation framework that systematically evaluates the accuracy of an assembly
    in a reference-independent manner using rigorous statistical methods.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize ALE2Wigglw
        :param camel: Camel instance
        :return: None
        """
        super().__init__('ALE2Wiggle', '2022.05.03', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        ale_input = self._folder / Path(str(self._tool_inputs['ALE'][0])).name
        ale_input.symlink_to(self._tool_inputs['ALE'][0].path)

        self._build_command(ale_input)
        self._execute_command()
        self._set_output()

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'ALE' not in self._tool_inputs:
            raise InvalidInputSpecificationError('ALE file is required')
        super()._check_input()

    def _build_command(self, ale_input) -> None:
        """
        Builds the command to run ALE2Wiggle.
        :return: None
        """
        self._command.command = ' '.join([self._tool_command,
                                          f'{ale_input}'])

    def _check_command_output(self) -> None:
        """
        Checks command output.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Command execution failed (Exit code: {self._command.returncode})")

    def _set_output(self) -> None:
        """
        Collects the tool outputs.
        """
        basename = Path(str(self._tool_inputs['ALE'][0])).name
        tsv_outputs = [f'{basename}-{stat}.wig' for stat in ['depth', 'kmer', 'insert', 'place']]
        self._tool_outputs['TSV'] = [ToolIOFile(self.folder / f'{f}') for f in tsv_outputs]
