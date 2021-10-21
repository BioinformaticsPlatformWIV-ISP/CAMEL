
import logging
import re
from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class PointFinder(Tool):
    """
    PointFinder performs detection of antimicrobial resistance associated with chromosomal point mutations in bacterial
    pathogens.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('PointFinder', '20190227', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError("FASTA input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._informs['last_update'] = self.__get_date_last_update()
        blastn_path = self.__determine_blastn_path()
        self.__build_command(blastn_path)
        self._execute_command()
        self.__set_output()
        self._informs['database'] = self._parameters['database'].value

    def __build_command(self, blastn_path: Path) -> None:
        """
        Builds the command line call.
        :param blastn_path: Path to the blastn binary
        :return: None
        """
        self._command.command = ' '.join([
            'source $VENV;',
            self._tool_command,
            '--inputfiles {}'.format(self._tool_inputs['FASTA'][0].path),
            '--method blastn',
            '--databasePath $POINTFINDER_DB',
            '--method_path {}'.format(blastn_path),
            '--out_path {}'.format(self._folder),
            ' '.join(self._build_options())
            ])

    def __determine_blastn_path(self) -> Path:
        """
        Determines the blastn path based on the dependencies.
        :return: blastn path
        """
        logging.debug('Retrieving blastn path')
        self._command.command = 'echo $PATH'
        self._execute_command()
        for path in self._command.stdout.split(':'):
            m = re.match('.*/blast/.*/bin', path)
            if m:
                return Path(path, 'blastn')
        raise ValueError("Cannot determine blastn path")

    def __get_date_last_update(self) -> str:
        """
        Returns the date of the last update.
        :return: Date of the last update.
        """
        logging.debug("Retrieving date of the last update.")
        self._command.command = 'echo $POINTFINDER_DB'
        self._execute_command()
        db_path = Path(self._command.stdout.strip())
        last_update_file = db_path / '.last_update.txt'
        if not last_update_file.is_file():
            raise FileNotFoundError('Cannot retrieve last update file: {}'.format(last_update_file))
        with last_update_file.open() as handle:
            return handle.readline().strip()

    def __set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        output_dir = Path(self._folder)
        self._tool_outputs['TSV'] = [ToolIOFile(next(
            f for f in output_dir.iterdir() if f.name.endswith('_results.tsv')))]
        self._tool_outputs['TXT'] = [ToolIOFile(next(
            f for f in output_dir.iterdir() if f.name.endswith('_HTMLtable.txt')))]

    def _check_command_output(self) -> None:
        """
        Checks the command output to see if the command executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Error executing command!")
