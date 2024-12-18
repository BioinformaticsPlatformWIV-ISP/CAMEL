from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class MashScreen(Tool):
    """
    Mash screen determines whether query sequences are within a larger mixture of sequences.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        :return: None
        """
        super().__init__('mash screen', '2.3', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'DB' not in self._tool_inputs:
            raise InvalidInputSpecificationError('Database input is required')
        if 'FASTQ' not in self._tool_inputs and 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('FASTQ or FASTA input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        path_tsv_out = self.folder / 'mash_screen.tsv'
        self._command.command = ' '.join([
            self._tool_command,
            str(self._tool_inputs['DB'][0].path),
            ' '.join(str(x.path) for x in self._tool_inputs['FASTQ']) if 'FASTQ' in self._tool_inputs else str(self._tool_inputs['FASTA'][0].path),
            *self._build_options(),
            f'> {path_tsv_out}'
        ])
        self._execute_command()
        self._tool_outputs['TSV'] = [ToolIOFile(path_tsv_out)]
        self._informs['cols'] = ['identity', 'hashes', 'm_mult', 'pval', 'q_id', 't']

    def _check_command_output(self) -> None:
        """
        Checks if the tool executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f'Error executing {self.name}: {self.stderr}')
