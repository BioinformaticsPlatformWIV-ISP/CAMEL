from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Unicycler(Tool):
    """
    Unicycler is an assembly pipeline for bacterial genomes.
    """

    OUTPUT_NAME = 'assembly.fasta'

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the Unicycler tool.
        :param camel: CAMEL instance
        """
        super().__init__('Unicycler', '0.5.0', camel)
        self._input_string = ''

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid.
        :return: None
        """
        if not any(key in ('FASTA', 'FASTQ_PE', 'FASTQ_SE') for key in self._tool_inputs):
            raise InvalidInputSpecificationError('FASTA or FASTQ_SE or FASTQ_PE input is required')
        super()._check_input()

    def _build_command(self) -> None:
        """
        Builds the command to run Unicycler.
        :return: None
        """
        if 'FASTQ_PE' in self._tool_inputs:
            self._input_string += f'--short1 {self._tool_inputs["FASTQ_PE"][0].path} ' \
                                  f'--short2 {self._tool_inputs["FASTQ_PE"][1].path} '
        if 'FASTQ_SE' in self._tool_inputs:
            self._input_string += f'--long {self._tool_inputs["FASTQ_SE"][0].path} '
        else:
            raise InvalidInputSpecificationError('FASTQ_SE input is required')
        self._command.command = ' '.join([
            self._tool_command, self._input_string, *self._build_options()])

    def _check_command_output(self) -> None:
        """
        Checks command output.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Command execution failed (Exit code: {self._command.returncode})")

    def _set_output(self, fasta_output: Path) -> None:
        """
        Collects the tool output.
        """
        dir_out = self.folder / self._parameters['output_dir'].value
        self._tool_outputs['FASTA'] = [ToolIOFile(dir_out / f'{fasta_output}')]

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        fasta_output = Path(Unicycler.OUTPUT_NAME)
        self._build_command()
        self._execute_command()
        self._set_output(fasta_output)
