from camel.app.camel import Camel
from camel.app.error.invalidparametererror import InvalidParameterError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.samtools.samtools import Samtools


class SamtoolsView(Samtools):
    """
    SAM <-> BAM Conversion
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super().__init__('samtools view', '1.9', camel)
        self.__input_key = None

    def _check_input(self) -> None:
        """
        Checks the input.
        :return: None
        """
        if 'SAM' in self._tool_inputs:
            self.__input_key = 'SAM'
        elif 'BAM' in self._tool_inputs:
            self.__input_key = 'BAM'
        else:
            raise ValueError("No input file found")
        super(Samtools, self)._check_input()

    def _check_parameters(self) -> None:
        """
        Checks the tool parameters.
        :return: None
        """
        if self._parameters['output_format'].value.upper() not in ('SAM', 'BAM'):
            raise InvalidParameterError("Invalid output format (BAM/SAM supported)")
        super(SamtoolsView, self)._check_parameters()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()
        self._check_stderr()

    def __build_command(self, pipe_in: bool = False, pipe_out: bool = False) -> None:
        """
        Builds the command for this tool.
        :return: None
        """
        excluded_parameters = ['regions']

        # Do not specify output filename when piping output
        if pipe_out:
            excluded_parameters.append('output_filename')

        # Construct command
        command_parts = [
            self._tool_command,
            ' '.join(self._build_options(excluded_parameters))
        ]

        # Add input file
        if not pipe_in:
            command_parts.append(str(self._tool_inputs[self.__input_key][0].path))

        # Add regions (when specified)
        if 'regions' in self._parameters:
            command_parts.append(f'"{self._parameters["regions"].value}"')

        self._command.command = ' '.join(command_parts)

    def __set_output(self) -> None:
        """
        Sets the tool output.
        :return: None
        """
        output_path = self.folder / self._parameters['output_filename'].value
        if not output_path.is_file:
            raise ToolExecutionError(f"Expected {self._name} output not generated")
        output_key = self._parameters['output_format'].value.upper()
        self._tool_outputs[output_key] = [ToolIOFile(output_path)]

    def _check_stderr(self) -> None:
        """
        Validates the stderr.
        :return: None
        """
        if 'only works for indexed' in self.stderr:
            raise ToolExecutionError("Can only extract regions from indexed BAM files")
        super(SamtoolsView, self)._check_stderr()

    def _before_pipe(self, dir_, pipe_in: bool, pipe_out) -> None:
        """
        Prepares the command that will be piped.
        :param dir_: Running directory
        :param pipe_in: True if tool receives piped input
        :param pipe_out: True if tool generates piped output
        :return: None
        """
        self.__build_command(pipe_in, pipe_out)

    def _after_pipe(self, stderr: str, is_last_in_pipe: bool) -> None:
        """
        Performs the required steps after executing the tool as part of a pipe.
        :param stderr: Stderr for this command in the pipe
        :param is_last_in_pipe: Boolean to indicate if this is the last step in the pipe
        :return: None
        """
        if is_last_in_pipe:
            self.__set_output()
