import json
from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class SeqSero2(Tool):
    """
    Salmonella serotype prediction from genome sequencing data.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize tool.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('SeqSero2', '1.2.1', camel)

    def _execute_tool(self) -> None:
        """
        Execute the tool.
        :return: None
        """
        self._informs['_tag'] = self._tool_inputs['MODE'][0].value
        self.build_command(self._informs['_tag'])
        self._execute_command()
        self.__set_output()
        input_folder = self._tool_inputs['DIR'][0].path
        self.__add_informs(input_folder)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        super(SeqSero2, self)._check_input()
        if not self._parameters.get('mode') \
                or self._parameters['mode'].value not in ('Kmer', 'Allele', 'Kmerread'):
            raise InvalidInputSpecificationError("A Seqsero2 processing mode must be passed to the tool, "
                                                 "choose from Kmer, Allele, or Kmerread")
        if self._parameters['mode'].value == 'Kmer':
            if 'FASTA' not in self._tool_inputs:
                raise InvalidInputSpecificationError("FASTA input is required in Kmer mode")
        else:
            if sum(x in self._tool_inputs for x in ('FASTQ', 'FASTQ_PE')) != 1:  # not exactly one
                raise InvalidInputSpecificationError(f"Exactly one FASTQ input is required in "
                                                     f"{self._parameters['mode'].value} mode")

    def __set_output(self) -> None:
        """
        Collects the tool output files.
        :return: None
        """
        self._tool_outputs['TXT'] = [ToolIOFile(self.folder / 'SeqSero_result.txt')]

    def build_command(self, mode: str) -> None:
        """
        Concatenates required parameters and options to build the command.
        :param mode: seqsero2 execution mode; 'Allele' or 'Kmer'
        :return: None
        """
        command_parts = [self._tool_command, '-d', str(self.folder), " ".join(self._build_options())]
        if mode == 'Kmer':
            command_parts.extend(['-t 4 -m k', '-i', str(self._tool_inputs['FASTA'][0])])
        else:
            if mode == 'Allele':
                command_parts.append('-m a')
            else:  # if self._tool_inputs['MODE'][0].value == 'Kmerread':
                command_parts.append('-m k')

            if 'FASTQ' in self._tool_inputs:
                command_parts.extend(['-t 3', '-i', str(self._tool_inputs['FASTQ'][0].path)])
            else:  # if 'FASTQ_PE' in self._tool_inputs:
                command_parts.extend(['-t 2', '-i',
                                      str(self._tool_inputs['FASTQ_PE'][0].path),
                                      str(self._tool_inputs['FASTQ_PE'][1].path)])
        self._command.command = ' '.join(command_parts)

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
        self._informs['db_path'] = input_folder
