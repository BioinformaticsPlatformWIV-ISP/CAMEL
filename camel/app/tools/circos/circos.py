from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Circos(Tool):
    """
    Circos is a software package for visualizing data and information. It visualizes data in a circular layout — this
    makes Circos ideal for exploring relationships between objects or positions.

    INPUT
        TXT: Text file with the configuration for Circos.
    OUTPUT
        PNG: Circos image rendered as PNG
        SVG: Circos image rendered as SVG
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL
        """
        super().__init__('Circos', '0.69-6', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is correct.
        :return: None
        """
        if 'TXT' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Circos configuration file ('TXT') is required.")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __build_command(self) -> None:
        """
        Builds the command line call.
        :return: None
        """
        self._command.command = ' '.join([self._tool_command, '-conf', str(self._tool_inputs['TXT'][0].path)])

    def __set_output(self) -> None:
        """
        Sets the tool output.
        :return: None
        """
        png_path = self.folder / 'circos.png'
        if not png_path.is_file():
            raise ToolExecutionError('No PNG output file generated')
        self._tool_outputs['PNG'] = [ToolIOFile(png_path)]
        svg_path = self.folder / 'circos.svg'
        if not svg_path.is_file():
            raise ToolExecutionError('No SVG output file generated')
        self._tool_outputs['SVG'] = [ToolIOFile(svg_path)]

    def _check_command_output(self) -> None:
        """
        Checks if the tool ran successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Error executing '{self.name}': {self._command.stderr}")
