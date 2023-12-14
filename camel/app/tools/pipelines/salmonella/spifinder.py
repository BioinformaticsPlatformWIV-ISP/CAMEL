import json
from pathlib import Path

from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.camel import Camel
from camel.app.tools.tool import Tool
from camel.app.io.tooliofile import ToolIOFile


class SPIFinder(Tool):

    """
    BLAST-based methodology tool for detection of SPI (Salmonella Pathogenicity Islands).
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize tool.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('SPIFinder', '0.1', camel)

    def _execute_tool(self) -> None:
        """
        Executes the tool
        :return: None
        """
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
        super(SPIFinder, self)._check_input()
        if not any(key in self._tool_inputs for key in ('FASTQ', 'FASTQ_PE', 'FASTA')):
            raise InvalidInputSpecificationError("FASTQ, FASTQ_PE or FASTA input is required")
        elif 'DIR' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Missing database!")

    def __set_output(self) -> None:
        """
        Sets the name of the output files
        :return: None
        """

        self._tool_outputs['JSON'] = [ToolIOFile(self.folder / 'data.json')]

    def __build_command(self) -> None:
        """
        Concatenates required parameters and options to build the command
        :return: None
        """

        if 'FASTA' in self._tool_inputs:
            self._informs['_tag'] = 'FASTA'
            self._command.command = ' '.join([
                self._tool_command, '-i',
                str(self._tool_inputs['FASTA'][0]), "-p ", str(self._tool_inputs['DIR'][0].path)
            ])
        elif 'FASTQ' in self._tool_inputs:
            self._informs['_tag'] = 'FASTQ'
            self._command.command = ' '.join([
                self._tool_command, '-i',
                str(self._tool_inputs['FASTQ'][0].path), "-p ", str(self._tool_inputs['DIR'][0].path)
            ])
        else:
            self._informs['_tag'] = 'FASTQ'
            self._command.command = ' '.join([
                self._tool_command, '-i',
                str(self._tool_inputs['FASTQ_PE'][0].path), str(self._tool_inputs['FASTQ_PE'][1].path),
                "-p ", str(self._tool_inputs['DIR'][0].path)
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
