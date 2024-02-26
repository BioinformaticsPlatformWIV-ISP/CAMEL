import json
from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Genotyphi(Tool):
    """
    Genotyphi scheme wrapped in a tool called Mykrobe (therefore the name of this class is poorly chosen, but
    renaming it would have too many implications in the bigsdb code & mongodb)
    The GenoTyphi genotyping scheme divides the Salmonella Typhi population into 4 major lineages, and >75 different
    clades and subclades.
    In addition, the tool also looks for mutations in some regions that confer resistance to
    several antibiotics families.
    """
    def __init__(self, camel: Camel) -> None:
        """
        Initializes Mykrobe.
        :param camel: Camel instance
        """
        super().__init__('mykrobe', 'v0.10.0', camel)

    def _execute_tool(self) -> None:
        """
        Runs Mykrobe with the Genotyphi scheme.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()
        input_folder = self._tool_inputs['DIR'][0].path
        self.__add_informs(input_folder)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        super(Genotyphi, self)._check_input()
        if 'DIR' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Database path needs to be specified")
        if ('FASTQ' in self._tool_inputs) == ('FASTQ_PE' in self._tool_inputs):  # not exactly one:
            raise InvalidInputSpecificationError("FASTQ or FASTQ_PE input is required")

    def __build_command(self) -> None:
        """
        Builds the command line call to execute Mykrobe with the Genotyphi scheme.
        :return: None
        """
        input_str = str(self._tool_inputs['FASTQ'][0].path) if 'FASTQ' in self._tool_inputs else \
            f"{str(self._tool_inputs['FASTQ_PE'][0].path)} {str(self._tool_inputs['FASTQ_PE'][1].path)}"
        self._command.command = ' '.join([
            self._tool_command, 'predict',
            '--sample Sample',
            '--species typhi',
            '--format csv',
            f'--seq {input_str}',
            *self._build_options()
        ])

    def _check_command_output(self) -> None:
        """
        Checks if the command output is valid.
        :return: None
        """
        if not self._command.returncode == 0:
            raise ToolExecutionError(f"Error executing {self.name}: {self._command.stderr}")

    def __set_output(self) -> None:
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs['CSV'] = [ToolIOFile(self.folder / str(self._parameters['output_filename'].value))]

    def __add_informs(self, input_folder: Path) -> None:
        """
        Adds the informs by parsing the JSON file containing the metadata in the database directory.
        :param input_folder: Input database directory
        :return: None
        """
        path_metadata = input_folder / 'db_update_info.json'
        if not path_metadata.is_file():
            raise FileNotFoundError(f'Database metadata not found: {path_metadata}')
        with path_metadata.open('r') as handle:
            metadata = json.load(handle)
        self._informs.update(metadata)
        