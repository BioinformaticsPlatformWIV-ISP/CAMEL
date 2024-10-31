from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.tools.toolpipeable import ToolPipeable


class Star(ToolPipeable):
    """
    Super class for STAR indexing of reference genomes and alignment of spliced transcripts
    """

    def __init__(self, tool_name: str, version: str, camel: Camel) -> None:
        """
        Initializes a STAR tool.
        :param tool_name: Tool name
        :param version: Tool version
        :param camel: Camel instance
        :return: None
        """
        super().__init__(tool_name, version, camel)
        self._required_inputs = []
        self._input_string = ""
        self._output_string = ""

    def _check_input(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        for key in self._required_inputs:
            if key not in self._tool_inputs:
                raise InvalidInputSpecificationError(f"Input '{key}' is required")
        super()._check_input()

    def _set_input(self) -> None:
        """
        Set the input specification
        :return: None
        """
        raise NotImplementedError("Method should be implemented by subclass.")

    def _build_command(self) -> None:
        """
        Builds the command to run STAR
        :return: None
        """
        self._command.command = ' '.join([self._tool_command, self._input_string,
                                          *self._build_options(excluded_parameters=['filename_output']),
                                          self._output_string])

    def _execute_tool(self) -> None:
        """
        Executes STAR
        :return: None
        """
        self._set_input()
        self._set_output()
        self._build_command()
        self._execute_command()

    def _set_output(self) -> None:
        """
        Set the output specification
        :return: None
        """
        raise NotImplementedError("Method should be implemented by subclass.")

    def _check_command_output(self) -> None:
        """
        Parse stderr message of STAR to check whether it runs successfully
        :return: None
        """
        if any(err in self.stderr.lower() for err in ('error', 'fail')):
            raise ToolExecutionError(f'Command failed: {self._command.command}\n stderr: {self.stderr}')