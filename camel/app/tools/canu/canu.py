from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Canu(Tool):
    """
    Canu is a fork of the Celera Assembler, designed for high-noise single-molecule sequencing (such as the PacBio RS
    II/Sequel or Oxford Nanopore MinION).
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Canu', '2.2 commit 7fb66bbff', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTQ' not in self._tool_inputs:
            raise InvalidInputSpecificationError("FASTQ input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            f"-nanopore-raw {self._tool_inputs['FASTQ'][0].path}",
            *self._build_options(['genome_size', 'threads', 'minimum_input_coverage']),
            *self._build_options(['output_directory', 'output_prefix'], '='),
            'useGrid=False',
            'stopOnLowCoverage=1'
        ])
        self._execute_command()
        self.__set_output()

    def _check_command_output(self) -> None:
        """
        Checks if the command executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Error executing '{self.name}': {self.stdout}")

    def __set_output(self) -> None:
        """
        Sets the output of the tool.
        :return: None
        """
        dir_out = self.folder / self._parameters['output_directory'].value
        self._tool_outputs['FASTA'] = [
            ToolIOFile(dir_out / f"{self._parameters['output_prefix'].value}.contigs.fasta")]
