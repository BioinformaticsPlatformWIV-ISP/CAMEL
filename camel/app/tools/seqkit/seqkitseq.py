from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class SeqkitSeq(Tool):

    """
    Seqkit seq performs common transformations of FASTA / FASTQ files.
    """

    INPUT_KEYS = ('FASTQ', 'FASTA')

    def __init__(self, camel: Camel) -> None:
        """
        Initialize seqkit seq.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Seqkit seq', '2.3.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if not any(x in self._tool_inputs for x in SeqkitSeq.INPUT_KEYS):
            raise InvalidInputSpecificationError('{} input is required.'.format(' or '.join(SeqkitSeq.INPUT_KEYS)))
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        output_path = self._folder / self._parameters['output_filename'].value
        self.__build_command()
        self._execute_command()
        self.__set_output(output_path)

    def __build_command(self) -> None:
        """
        Builds the command line call.
        :return: None
        """
        input_key = 'FASTQ' if 'FASTQ' in self._tool_inputs else 'FASTA'
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join([str(f.path) for f in self._tool_inputs[input_key]]),
            *self._build_options()
        ])

    def __set_output(self, output_path: Path) -> None:
        """
        Sets the tool output.
        :param output_path: Output path
        :return: None
        """
        output_key = 'FASTA' if output_path.name.lower().endswith('.fasta') else 'FASTQ'
        self._tool_outputs[output_key] = [ToolIOFile(output_path)]

    def _check_command_output(self) -> None:
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Command execution failed (Exit code: {})".format(self._command.returncode))
