import logging

import os
import re

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

    def __init__(self, camel: Camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('PointFinder', '3.0', camel)

    def _check_input(self):
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

    def __build_command(self, blastn_path: str) -> None:
        """
        Builds the command line call.
        :param blastn_path: Path to the blastn binary
        :return: None
        """
        self._command.command = ' '.join([
            'source $VENV;',
            self._tool_command,
            '--inputfiles {}'.format(self._tool_inputs['FASTA'][0].path),
            '--databasePath $POINTFINDER_DB',
            '--blastPath {}'.format(blastn_path),
            '--out_path {}'.format(self._folder),
            ' '.join(self._build_options())
            ])

    def __determine_blastn_path(self) -> str:
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
                return os.path.join(path, 'blastn')
        raise ValueError("Cannot determine blastn path")

    def __get_date_last_update(self) -> str:
        """
        Returns the date of the last update.
        :return: Date of the last update.
        """
        logging.debug("Retrieving date of the last update.")
        self._command.command = 'echo $POINTFINDER_DB'
        self._execute_command()
        db_path = self._command.stdout.strip()
        last_update_file = os.path.join(db_path, '.last_update.txt')
        if not os.path.isfile(last_update_file):
            raise FileNotFoundError('Cannot retrieve last update file: {}'.format(last_update_file))
        with open(last_update_file) as handle:
            return handle.readline().strip()

    def __set_output(self) -> None:
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['TSV'] = [ToolIOFile(os.path.join(self._folder, 'PointFinder_results.txt'))]
        self._tool_outputs['TXT'] = [ToolIOFile(os.path.join(self._folder, 'PointFinder_table.txt'))]

    def _check_command_output(self):
        """
        Checks the command output to see if the command executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Error executing command!")
