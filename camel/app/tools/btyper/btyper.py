import os
from pathlib import Path
from camel.app.camel import Camel
from camel.app.tools.tool import Tool
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError

class BTyper(Tool):


    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('BTyper', '3.2.0', camel)

    def _check_input(self) -> None:
        """
        Checks whether the required input files are specified.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('No FASTA input found')
        super()._check_input()

        my_path = self._folder / Path(str(self._tool_inputs['FASTA'][0])).name
        os.symlink(str(self._tool_inputs["FASTA"][0]), str(my_path))
        self._tool_inputs['FASTA'][0] = my_path

    def _build_command(self) -> None:
        """
        Build the command to run tool
        :return: None
        """
        self._command.command = f'{self._tool_command} ' \
                                f'--input {self._tool_inputs["FASTA"][0]} ' \
                                f'{" ".join(self._build_options())}'

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._build_command()
        self._execute_command()
        self._tool_outputs['TSV'] = [ToolIOFile(self._parameters['output_dir'].value / 'btyper3_final_results' / 'contigs_final_results.txt')]