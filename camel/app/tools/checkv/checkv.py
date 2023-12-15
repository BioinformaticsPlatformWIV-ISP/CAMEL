from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class CheckV(Tool):
    """
    Assessing the quality of metagenome-assembled viral genomes
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('CheckV', '1.0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('FASTA input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        dir_out = self.folder / 'out'
        self._command.command = ' '.join([
            self._tool_command,
            'end_to_end',
            str(self._tool_inputs['FASTA'][0].path),
            str(dir_out),
            *self._build_options()
        ])
        self._execute_command()
        self.__set_output(dir_out)

    def __set_output(self, dir_out: Path) -> None:
        """
        Sets the tool output for this tool.
        :param dir_out: Output directory
        :return: None
        """
        for filename in [Path(x) for x in (
                'complete_genomes.tsv', 'completeness.tsv', 'contamination.tsv', 'quality_summary.tsv')]:
            self._tool_outputs[f'TSV_{filename.stem}'] = [ToolIOFile(dir_out / filename)]

    def _check_command_output(self) -> None:
        """
        Checks if the tool was executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Error executing {self.name}: {self._command.stderr}")
