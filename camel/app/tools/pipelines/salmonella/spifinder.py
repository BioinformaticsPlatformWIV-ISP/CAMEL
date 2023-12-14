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
        db_dir = self._tool_inputs['DIR'][0].path
        self.__add_informs(db_dir)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        super(SPIFinder, self)._check_input()
        # check if exactly one of the three possible inputs is provided
        input_keys = [key for key in ('FASTQ', 'FASTQ_PE', 'FASTA') if key in self._tool_inputs]
        if len([key for key in ('FASTQ', 'FASTQ_PE', 'FASTA') if key in self._tool_inputs]) != 1:
            raise InvalidInputSpecificationError("Exactly one of FASTQ, FASTQ_PE or FASTA input is required.")
        else:
            self._input_key = input_keys[0]
        if 'DIR' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Database input is required (DIR).")

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
        self._informs['_tag'] = 'FASTQ' if self._input_key != 'FASTA' else 'FASTA'
        if self._input_key == 'FASTQ_PE':
            inputs_str = ' '.join([
                str(self._tool_inputs[self._input_key][0].path),
                str(self._tool_inputs[self._input_key][1].path)
            ])
        else:
            inputs_str = str(self._tool_inputs[self._input_key][0])

        self._command.command = ' '.join([
            self._tool_command, '-i', inputs_str, "-p", str(self._tool_inputs['DIR'][0].path)
        ])

    def _check_command_output(self) -> None:
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if 'error' in self.stderr.lower() or self._command.returncode != 0:
            raise ToolExecutionError(f"Error executing {self.name}: {self._command.stderr.strip()}")

    def __add_informs(self, db_dir: Path) -> None:
        """
        Adds the informs by parsing the JSON file containing the metadata in the database directory.
        :param db_dir: Input database directory
        :return: None
        """
        self._informs['_tag'] = 'FASTQ' if not 'FASTA' in self._tool_inputs else 'FASTA'
        db_metadata_file = db_dir / 'db_update_info.json'
        if not db_metadata_file.is_file():
            raise FileNotFoundError(f'Database metadata file not found: {db_metadata_file}')
        with db_metadata_file.open('r') as handle:
            metadata = json.load(handle)
        self._informs.update(metadata)
