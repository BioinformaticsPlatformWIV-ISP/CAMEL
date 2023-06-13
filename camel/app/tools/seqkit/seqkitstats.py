import pandas as pd

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class SeqkitStats(Tool):
    """
    Seqkit stats reports simple statistics of FASTA/Q files.
    """

    INPUT_KEYS = ('FASTQ', 'FASTA')

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Seqkit stats', '2.3.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if not any(x in self._tool_inputs for x in SeqkitStats.INPUT_KEYS):
            raise InvalidInputSpecificationError('{} input is required.'.format(' or '.join(SeqkitStats.INPUT_KEYS)))
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        tsv_out = self.folder / self._parameters['output_filename'].value
        self.__build_command()
        self._execute_command()
        data = pd.read_table(tsv_out)
        for k, v in next(iter(data.to_dict('records'))).items():
            self._informs[k] = v
        self._tool_outputs['TSV'] = [ToolIOFile(tsv_out)]

    def __build_command(self) -> None:
        """
        Builds the command line call.
        :return: None
        """
        input_key = 'FASTQ' if 'FASTQ' in self._tool_inputs else 'FASTA'
        self._command.command = ' '.join([
            self._tool_command,
            '--all', '--tabular',
            ' '.join([str(f.path) for f in self._tool_inputs[input_key]]),
            *self._build_options()
        ])

    def _check_command_output(self) -> None:
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Command execution failed (Exit code: {})".format(self._command.returncode))
