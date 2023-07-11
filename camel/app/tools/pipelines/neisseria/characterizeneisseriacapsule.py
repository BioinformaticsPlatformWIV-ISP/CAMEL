import os
import shutil
import tempfile
from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class CharacterizeNeisseriaCapsule(Tool):
    """
    characterize_neisseria_capsule is a tool implementing a WGS-based method for N. meningitidis
    serogroup predictions by identifying capsule genes and genetic variations that might impact their expression.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the characterize_neisseria_capsule tool.
        :param camel: CAMEL instance
        """
        super().__init__('characterize_neisseria_capsule', 'a75a009', camel)

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid:
        - FASTA_dir is the only required input
        :return: None
        """
        if 'FASTA_dir' not in self._tool_inputs:
            raise InvalidInputSpecificationError('Required input directory containing the fasta files is missing')
        super()._check_input()

    def __build_input_string(self, dir_path: Path, dir_out: Path) -> str:
        """
        Creates the string with the input files
        :return: String with the input parameters
        """
        inputs = [f"-o {dir_out}",
                  f"-t {self._parameters['threads'].value}",
                  f"-d {dir_path}"]
        return ' '.join(inputs)

    def __build_command(self, dir_path: Path, dir_out: Path) -> None:
        """
        Concatenates required parameters and options to build the command.
        :return: None
        """
        input_string = self.__build_input_string(dir_path, dir_out)
        self._command.command = ' '.join([self._tool_command, input_string])

    def _check_command_output(self) -> None:
        """
        Checks if the command executed successfully.
        :return: None
        """
        if not self._command.returncode == 0:
            raise ToolExecutionError(f'Error executing {self.name}: {self.stderr}')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        temp_input_dir = tempfile.TemporaryDirectory()
        try:
            for filename in os.scandir(self._tool_inputs['FASTA_dir'][0].path):
                if os.path.splitext(filename)[-1] == '.fasta':
                    # Temporary copy to use a single FASTA file as input
                    shutil.copy(filename.path, Path(temp_input_dir.name))

                    # Execute the tool
                    sample_name = os.path.splitext(os.path.split(filename)[-1])[0]
                    dir_out = Path(f"{self._parameters['output_directory'].value}/{sample_name}")
                    self.__build_command(temp_input_dir.name, dir_out)
                    self._execute_command()

                    # Remove temporary fasta file
                    filepath = os.path.join(temp_input_dir.name, os.path.basename(filename))
                    os.remove(filepath)
                try:
                    self._tool_outputs['TSV'] = [ToolIOFile(next((dir_out / 'serogroup').glob('serogroup_predictions_'
                                                                                              '*.tab')))]
                except StopIteration:
                    raise ToolExecutionError(f"TSV file not found in output folder: {dir_out / 'serogroup'}")
        finally:
            temp_input_dir.cleanup()