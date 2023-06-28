from pathlib import Path
from typing import List

from camel.app.camel import Camel
from camel.app.components.files.fastautils import FastaUtils
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Polypolish(Tool):

    """
    Polishing assembly with short reads using polypolish.

    INPUT:
    - FASTA file of assembly
    - SAM file containing read mapping information

    OUTPUT:
    - FASTA file of polished assembly
    """
    OUTPUT_NAME = 'polished.fasta'

    def __init__(self, camel: Camel) -> None:
        """
        Initializes Polypolish.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Polypolish', '0.5.0', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        fasta_input = self._tool_inputs['FASTA'][0].path
        sam_input = [sam.path for sam in self._tool_inputs['SAM']]
        fasta_output = Path(Polypolish.OUTPUT_NAME)
        self._build_command(fasta_input, sam_input, fasta_output)
        self._execute_command()
        self._set_output(fasta_output)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('FASTA reference is required')
        if 'SAM' not in self._tool_inputs:
            raise InvalidInputSpecificationError('SAM alignment file is required')

        if len(self._tool_inputs['SAM']) > 2:
            raise InvalidInputSpecificationError('Please input at most two SAM alignment files')
        if not FastaUtils.is_indexed(self._tool_inputs['FASTA'][0].path):
            raise InvalidInputSpecificationError('FASTA reference needs to be indexed')
        super()._check_input()

    def _build_command(self, fasta_input: Path, sam_input: List[Path], fasta_output: Path) -> None:
        """
        Builds the command to run Polypolish.
        :param fasta_input: Assembly to polish
        :param sam_input: list of SAM alignment files
        :param fasta_output: Polished assembly
        :return: None
        """
        self._command.command = ' '.join([self._tool_command,
                                          *self._build_options(),
                                          str(fasta_input),
                                          ' '.join(str(sam_file) for sam_file in sam_input),
                                          f'> {fasta_output}'])

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
        :param fasta_output: Path of the output fasta file
        :return: None
        """
        self._tool_outputs['FASTA'] = [ToolIOFile(self.folder / f'{fasta_output}')]
