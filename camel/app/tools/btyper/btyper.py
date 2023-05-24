from pathlib import Path
from camel.app.camel import Camel
from camel.app.tools.tool import Tool
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile


class BTyper(Tool):
    """
    In silico taxonomic classification of Bacillus cereus group isolates using assembled genomes.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('BTyper', '3.3.4', camel)

    def _check_input(self) -> None:
        """
        Checks whether the required input files are specified.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('No FASTA input found')
        super()._check_input()

    def _build_command(self, fasta_input: Path) -> None:
        """
        Build the command to run tool.
        :return: None
        """
        self._command.command = ' '.join([self._tool_command, f'--input {fasta_input}', *self._build_options()])

    def _check_command_output(self) -> None:
        """
        Checks command output.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Command execution failed (Exit code: {self._command.returncode})")

    def _set_output(self) -> None:
        """
        set the output file to check.
        """
        output_filename = f'btyper3_final_results/{self._tool_inputs["FASTA"][0].path.stem}_final_results.txt'
        self._tool_outputs['TSV'] = [ToolIOFile(Path(self._parameters['output_dir'].value) / Path(output_filename))]

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Symlink the input FASTA file
        fasta_input = self._folder / Path(str(self._tool_inputs['FASTA'][0])).name
        try:
            fasta_input.symlink_to(self._tool_inputs["FASTA"][0].path)
        except FileExistsError:
            pass

        # Building the command
        self._build_command(fasta_input)
        self._execute_command()

        # Collect output
        self._set_output()
