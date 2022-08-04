from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class ResFinder(Tool):
    """
    ResFinder identifies acquired antimicrobial resistance genes in total or partial sequenced isolates of bacteria.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool
        :param camel: CAMEL instance
        """
        super().__init__('ResFinder', '4.1.11', camel)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if any(key in self._tool_inputs for key in ('FASTA' or 'FASTQ_PE' or 'FASTQ_SE')):
            raise InvalidInputSpecificationError('FASTA or FASTQ input is required')
        super()._check_input()

    def _build_command(self) -> None:
        """
        Builds the command to run resfinder.
        :return: None
        """
        if 'FASTA' in self._tool_inputs:
            input_str = f'--inputfasta {self._tool_inputs["FASTA"][0].path}'
        elif 'FASTQ_SE' in self._tool_inputs:
            input_str = f'--inputfastq {self._tool_inputs["FASTQ"][0].path}'
        elif 'FASTQ_PE' in self._tool_inputs:
            input_str = f'--inputfastq {self._tool_inputs["FASTQ_PE"][0].path} ' \
             f'{str(self._tool_inputs["FASTQ_PE"][1].path)}'
        self._command.command = ' '.join([self._tool_command, input_str, ' '.join(self._build_options())])

    def _check_command_output(self) -> None:
        """
        Checks command output.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Command execution failed (Exit code: {self._command.returncode})")

    def _set_output(self) -> None:
        """
        Collects the tool output.
        """
        self._tool_outputs['TSV'] = [
            ToolIOFile(self._parameters['output_path'].value / Path('ResFinder_results_tab.txt'))]
        if self._parameters['point'] is True:
            self._tool_outputs['TSV'].append(
                ToolIOFile(self._parameters['output_path'].value / Path('PointFinder_results.txt'))
            )

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._build_command()
        self._execute_command()
        self._set_output()
