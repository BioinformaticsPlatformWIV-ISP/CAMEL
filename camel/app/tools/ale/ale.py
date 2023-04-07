from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class ALE(Tool):

    """
    ALE is the Assembly Likelihood Evaluation framework that systematically evaluates the accuracy of an assembly
    in a reference-independent manner using rigorous statistical methods.
    """

    ALE_OUTPUT = "ALE.ale"

    def __init__(self, camel: Camel) -> None:
        """
        Initialize Polypolish
        :param camel: Camel instance
        :return: None
        """
        super().__init__('ALE', '2022.05.03', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        fasta_input = Path(str(self._tool_inputs['FASTA'][0]))
        bam_input = Path(str(self._tool_inputs['SAM'][0]))
        ale_output = Path(ALE.ALE_OUTPUT)
        self._build_command(fasta_input, bam_input, ale_output)
        self._execute_command()
        self._set_output(ale_output)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('FASTA reference is required')
        if 'SAM' not in self._tool_inputs:
            raise InvalidInputSpecificationError('SAM alignment file is required')

        input_folder = self._tool_inputs['FASTA'][0].path.parent
        base_fasta_name = self._tool_inputs['FASTA'][0].path.name
        fasta_index_file = [f for f in input_folder.glob(f'{base_fasta_name}.fai')]
        if not (len(fasta_index_file) > 0):
            raise InvalidInputSpecificationError('FASTA reference needs to be indexed')
        super()._check_input()

    def _build_command(self, fasta_input: Path, sam_input: Path, ale_output: Path) -> None:
        """
        Builds the command to run ALE.
        :return: None
        """
        self._command.command = ' '.join([self._tool_command,
                                          *self._build_options(),
                                          f'{sam_input}',
                                          f'{fasta_input}',
                                          f'{ale_output}'])

    def _check_command_output(self) -> None:
        """
        Checks command output.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Command execution failed (Exit code: {self._command.returncode})")

    def _set_output(self, ale_output: Path) -> None:
        """
        Collects the tool output.
        """
        self._tool_outputs['ALE'] = [ToolIOFile(self.folder / f'{ale_output}')]
