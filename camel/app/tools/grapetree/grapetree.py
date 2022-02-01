from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class GrapeTree(Tool):
    """
    Package to create cgMLST-based phylogenies associated with EnteroBase.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('GrapeTree', '2.2', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Allele matrix (TSV) input is required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        output_path = self.folder / self._parameters['output_path'].value
        self._command.command = ' '.join([
            self._tool_command,
            '--profile', str(self._tool_inputs['TSV'][0].path)
        ] + self._build_options(excluded_parameters=['output_path']) + [f'> {output_path}'])
        self._execute_command()
        self._tool_outputs['NWK'] = [ToolIOFile(output_path)]

    def _check_command_output(self) -> None:
        """
        Checks if command executed successfully.
        :return: None
        """
        if not self._command.returncode == 0:
            raise ToolExecutionError(self._command.stderr)
