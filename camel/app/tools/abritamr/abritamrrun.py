import json
from pathlib import Path

from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.camel import Camel
from camel.app.tools.tool import Tool
from camel.app.io.tooliofile import ToolIOFile


class AbriTAMRRun(Tool):
    """
    AbritAMR: abriTAMR is an AMR gene detection pipeline that runs AMRFinderPlus on a single (or list ) of given
    isolates and collates the results into a table, separating genes identified into functionally relevant groups.
    This is the first part of the AbriTAMR pipeline (run).
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize tool.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('AbriTAMR run', '1.0.13', camel)

    def _execute_tool(self) -> None:
        """
        Executes the tool
        :return: None
        """
        self.__set_output()
        self.__build_command()
        self._execute_command()
        amrfinder_db_folder = self._tool_inputs['DIR_AMRF'][0].path
        self.__add_informs(amrfinder_db_folder)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        super(AbriTAMRRun, self)._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError("FASTA input is required")
        elif 'DIR_AMRF' not in self._tool_inputs:
            raise InvalidInputSpecificationError("databases path need to be specified")
        elif 'VAL_SPECIES' not in self._tool_inputs:
            raise InvalidInputSpecificationError(f"A species must be provided to run {self.name}")

    def __set_output(self) -> None:
        """
        Collects the output files of interest.
        :return: None
        """

        self._tool_outputs['TXT_MATCHES'] = [ToolIOFile(self.folder / 'summary_matches.txt')]
        self._tool_outputs['TXT_PARTIALS'] = [ToolIOFile(self.folder / 'summary_partials.txt')]

    def __build_command(self) -> None:
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        self._informs['_tag'] = 'RUN'
        self._command.command = ' '.join([
            self._tool_command, '--contigs',
            str(self._tool_inputs['FASTA'][0]), '--prefix', str(self.folder), '--species',
            str(self._tool_inputs['VAL_SPECIES'][0]),
            '--amrfinder_db', str(self._tool_inputs['DIR_AMRF'][0].path)
        ])

    def _check_command_output(self) -> None:
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if 'error' in self.stderr.lower():
            raise ToolExecutionError(f"Command execution failed (stderr: {self.stderr}).")
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Command execution failed (Exit code: {self._command.returncode})")

    def __add_informs(self, amrfinder_folder: Path) -> None:
        """
        add the update of the two database in the informs of the tool for further reporting
        amrfinder_folder: the path to the folder of amrfinder plus to be used with the tool abritAMR
        return: None
        """
        amrfinder_db_metadata = amrfinder_folder / 'db_update_info.json'
        if not amrfinder_db_metadata.is_file():
            raise FileNotFoundError(f'Database metadata not found: {amrfinder_db_metadata}')
        with amrfinder_db_metadata.open() as handle:
            metadata = json.load(handle)
            self._informs.update(metadata)
