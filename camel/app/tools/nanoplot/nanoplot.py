import pandas as pd

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class NanoPlot(Tool):
    """
    Plotting tool for long read sequencing data and alignments.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('NanoPlot', '1.36.2', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        """
        if 'FASTQ' not in self._tool_inputs:
            raise InvalidInputSpecificationError('FASTQ input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        """
        self._command.command = ' '.join([
            self._tool_command, '--fastq', str(self._tool_inputs['FASTQ'][0].path), *self._build_options()
        ])
        self._execute_command()
        self.__set_output()
        self.__set_informs()

    def __set_output(self) -> None:
        """
        Collects the tool output and stores it in tool outputs.
        """
        self._tool_outputs['TSV'] = [ToolIOFile(self.folder / 'NanoStats.txt')]
        self._tool_outputs['HTML'] = [ToolIOFile(self.folder / 'NanoPlot-report.html')]

    def _check_command_output(self) -> None:
        """
        Checks command output.
        :return: None
        """
        for line in self.stderr.splitlines():
            if 'skipping' in line:
                print(f"WARNING: {line}")
                self._command.returncode = 0
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Command execution failed (Exit code: {self._command.returncode})")

    def __set_informs(self) -> None:
        """
        Collects the informs.
        """
        data_summary = pd.read_table(self._tool_outputs['TSV'][0].path)
        for row in data_summary.to_dict('records'):
            self._informs[row['Metrics']] = row['dataset']
