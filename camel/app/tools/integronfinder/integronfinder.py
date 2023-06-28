import pandas as pd

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class IntegronFinder(Tool):
    """
    Finds integrons in DNA sequences.
    Reference: https://github.com/gem-pasteur/Integron_Finder
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('IntegronFinder', '2.0.2', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('FASTA input is required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            str(self._tool_inputs['FASTA'][0].path),
            *self._build_options()
        ])
        self._execute_command()
        dir_out = next(d for d in self.folder.glob('Results_*') if d.is_dir())
        path_tsv_out = next(dir_out.glob('*.integrons'))
        try:
            data_out = pd.read_table(path_tsv_out, comment='#')
            self._informs['nb_detected'] = len(data_out)
        except pd.errors.EmptyDataError as err:
            self._informs['nb_detected'] = 0
        self._tool_outputs['TSV'] = [ToolIOFile(path_tsv_out)]

    def _check_command_output(self) -> None:
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if not self._command.returncode == 0:
            raise ToolExecutionError(self._command.stderr)
