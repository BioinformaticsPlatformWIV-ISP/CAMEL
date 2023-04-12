from pathlib import Path

from camel.app.camel import Camel
from camel.app.components.files.fastautils import FastaUtils
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Sniffles(Tool):
    """
    Sniffles is a fast structural variant caller for long-read sequencing, it accurately detects SVs on germline,
    somatic and population-level for PacBio and Oxford Nanopore read data.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes Sniffles 2.0.7.
        :param camel: CAMEL instance
        :return: None
        """
        super().__init__('Sniffles', '2.0.7', camel)

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('FASTA file is required')
        if 'BAM' not in self._tool_inputs:
            raise InvalidInputSpecificationError('BAM alignment file is required')

        if not FastaUtils.is_indexed(self._tool_inputs['FASTA'][0].path, self._tool_inputs['FASTA'][0].path.parent):
            raise InvalidInputSpecificationError('FASTA reference needs to be indexed')
        super()._check_input()

    def _build_command(self, fasta_input: Path, bam_input: Path) -> None:
        """
        Builds the command to run sniffles.
        :return: None
        """
        self._command.command = ' '.join([self._tool_command,
                                          f'--input {bam_input}',
                                          f'--reference {fasta_input}',
                                          *self._build_options()])

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
        self._tool_outputs['VCF'] = [ToolIOFile(self.folder / Path('variants.vcf'))]

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        fasta_input = Path(str(self._tool_inputs['FASTA'][0]))
        bam_input = Path(str(self._tool_inputs['BAM'][0]))
        self._build_command(fasta_input, bam_input)
        self._execute_command()
        self._set_output()
