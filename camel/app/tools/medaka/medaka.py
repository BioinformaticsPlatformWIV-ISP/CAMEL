import abc
from pathlib import Path

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class Medaka(Tool, metaclass=abc.ABCMeta):

    """
    Super class for Medaka tools
    """

    def __init__(self, tool_name, version, camel):
        """
        Initialize a Medaka tool
        :param tool_name: Tool name
        :param version: Tool version
        :param camel: Camel instance
        :return: None
        """
        super().__init__(tool_name, version, camel)
        self._specific_parameters = []
        self._required_inputs = []
        self._input_string = ''
        self._output_string = ''
        self._option_string = ''
        self._output_type = ''

    def _execute_tool(self):
        """
        Run a Medaka function
        :return: None
        """
        self._set_input()
        self._set_output()
        self._set_specific_parameters()
        self._build_command()
        self._execute_command()
        self._set_informs()

    def _check_input(self):
        """
        Check input for a tool and prepare command line parameters for input
        :return: None
        """
        for input_file in self._required_inputs:
            if input_file not in self._tool_inputs:
                raise InvalidInputSpecificationError(
                    'Medaka {!r} required {!r} input is missing in _tool_inputs!'.format(self._name, input_file))

        super(Medaka, self)._check_input()

    def _set_input(self):
        """
        Set the input specification
        :return: None
        """
        pass

    def _set_output(self):
        """
        Set the output specification
        :return: None
        """
        self._output_string = "{}".format(self._parameters['output'])
        self._tool_outputs[self._output_type] = [ToolIOFile(Path(self._folder / self._parameters['output'].value))]

    def _set_specific_parameters(self):
        """
        Set specific parameters that need special handling when required
        :return: None
        """
        pass

    def _build_command(self):
        """
        Build the command to run tool
        :return: None
        """
        self._option_string += " ".join(self._build_options(excluded_parameters=['output']))
        self._command.command = " ".join([
            self._tool_command, self._option_string, self._input_string, self._output_string
        ])

    def _check_command_output(self) -> None:
        """
        Checks command output.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Command execution failed (Exit code: {self._command.returncode})")

    def _set_informs(self):
        """
        Set informs by analyze the output
        :return: None
        """
        pass
