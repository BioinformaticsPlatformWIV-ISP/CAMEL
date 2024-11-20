from os import listdir, path
from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class PharokkaMultiplotter(Tool):
    """
    Pharokka Multiplotter is a companion functionality of Pharokka to visualize genome annotations.
    """

    OUTPUT_DICT = {
        'PNG': '.png',
        'SVG': '.svg'
    }

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super().__init__('Pharokka', '1.7.3', camel)

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()

        # Collect the output
        path_out = Path(self.folder, self._parameters['outdir'].value)
        plots = listdir(path_out)
        extensions = PharokkaMultiplotter.OUTPUT_DICT.values()

        for file in plots:
            file_extension = path.splitext(file)[1]
            if file_extension not in extensions:
                raise ToolExecutionError(f'{file} not found in {path_out}')
            file_key = [key for key, ext in PharokkaMultiplotter.OUTPUT_DICT.items() if ext == file_extension][0]
            try:
                self._tool_outputs[file_key] = self._tool_outputs[file_key] + [ToolIOFile(Path(path_out, file))]
            except:
                self._tool_outputs[file_key] = [ToolIOFile(Path(path_out, file))]

    def _check_input(self) -> None:
        """
        Checks the tool input.
        :return: None
        """
        if 'GBK' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No genbank input found")
        if len(self._tool_inputs['GBK']) != 1:
            raise InvalidInputSpecificationError("Only one genbank file can be annotated at a time.")
        super()._check_input()

    def _check_command_output(self) -> None:
        """
        Checks the command output to see if the program ran correctly.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Error executing Pharokka Multiplotter: {self.stderr}")

    def __build_command(self) -> None:
        """
        Build the command line call.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            *self._build_options(),
            '-g',
            str(self._tool_inputs['GBK'][0].path)
        ])
