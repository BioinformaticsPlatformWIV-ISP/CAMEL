from camel.app.camel import Camel
from camel.app.components.files.fastautils import FastaUtils
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Polca(Tool):
    """
    Polishing assembly with short reads using polca.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes Polca.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('POLCA', '4.1.0', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._build_command()
        self._execute_command()
        self._set_output()

    def _check_input(self) -> None:
        """
        Checks whether the provided input files are valid
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('FASTA reference is required')
        if 'FASTQ_PE' not in self._tool_inputs:
            raise InvalidInputSpecificationError('FASTQ_PE files are required')

        if not FastaUtils.is_indexed(self._tool_inputs['FASTA'][0].path):
            raise InvalidInputSpecificationError('FASTA reference needs to be indexed')
        super()._check_input()

    def _build_command(self) -> None:
        """
        Builds the command to run polca.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            f"-a {self._tool_inputs['FASTA'][0]}",
            f"-r {self._tool_inputs['FASTQ_PE'][0]}",
            f"-r {self._tool_inputs['FASTQ_PE'][1]}",
            *self._build_options()])

    def _output_is_empty(self) -> bool:
        """
        Function to check if the VCF output is empty (no variants detected).
        :return: True if output is empty
        """
        return 'fasta.vcf: No such file or directory' in self._command.stderr

    def _check_command_output(self) -> None:
        """
        Checks command output.
        :return: False or None
        """
        if self._output_is_empty():
            return
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Command execution failed (Exit code: {self._command.returncode})")

    def _set_output(self) -> None:
        """
        Collects the tool output.
        """
        fasta_output = self.folder / f'{self._tool_inputs["FASTA"][0].path.name}.PolcaCorrected.fa'
        if self._output_is_empty():
            fasta_output.symlink_to(self._tool_inputs['FASTA'][0].path)
        self._tool_outputs['FASTA'] = [ToolIOFile(fasta_output)]
