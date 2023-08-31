from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Minimap2Mapping(Tool):
    """
    A versatile pairwise aligner for genomic and spliced nucleotide sequences.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool
        :param camel: CAMEL instance
        """
        super().__init__('Minimap2', '2.26', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        path_out = self.folder / self._parameters['output_filename'].value
        self.__build_command(path_out)
        self._execute_command()
        self.__set_output(path_out)

    def __build_command(self, path_out: Path) -> None:
        """
        Builds the command line call.
        :param path_out: Output filename
        :return: None
        """
        parts = [
            self._tool_command,
            '-ax map-ont',
            str(self._tool_inputs['FASTA'][0].path),
            str(self._tool_inputs['FASTQ'][0].path),
            ' '.join(self._build_options(excluded_parameters=['output_filename'])),
            f'> {path_out}'
        ]
        self._command.command = ' '.join(parts)

    def _check_command_output(self) -> None:
        """
        Checks if the tool executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Error executing {self.name}:\n{self.stderr}")

    def __set_output(self, path_out: Path) -> None:
        """
        Sets the output of this tool.
        :param path_out: Output filename
        :return: None
        """
        self._tool_outputs['SAM'] = [ToolIOFile(Path(path_out))]
