from pathlib import Path
import json

from camel.app.camel import Camel
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Sistr(Tool):

    """
    Serovar predictions from whole-genome sequence assemblies by determination of antigen
    gene and cgMLST gene alleles using BLAST.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize tool.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('SISTR', '1.1.1', camel)

    def _execute_tool(self):
        self.__set_output()
        self.__build_command()
        self._execute_command()
        input_folder = self._tool_inputs['DIR'][0].path
        self.__add_informs(input_folder)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        super(Sistr, self)._check_input()
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError("FASTA input is required")

    def __set_output(self) -> None:
        """
        Sets the name of the output files
        :return: None
        """

        self._tool_outputs['TSV'] = [ToolIOFile(self.folder / 'sistr_output.json')]

    def __build_command(self) -> None:
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command, '-f json --use-full-cgmlst-db --qc -v -o', str(self.folder / 'sistr_output.json'),
            str(self._tool_inputs['FASTA'][0].path), ' '.join(self._build_options())
            ])

    def _check_command_output(self) -> None:
        """
        Checks if the command was executed successfully.
        Here doesn't check for errors as some small warnings are displayed by the tool as error so not meaningful.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Command execution failed (Exit code: {self._command.returncode})")

    def __add_informs(self, input_folder: Path) -> None:
        """
        Adds the informs by parsing the JSON file containing the metadata in the database directory.
        :param input_folder: Input database directory
        :return: None
        """
        path_metadata = input_folder / 'db_update_info.json'
        if not path_metadata.is_file():
            raise FileNotFoundError(f'Database metadata not found: {path_metadata}')
        with path_metadata.open() as handle:
            metadata = json.load(handle)
        self._informs.update(metadata)
        self._informs['db_path'] = input_folder
