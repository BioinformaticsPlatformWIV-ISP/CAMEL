from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class NcbiHumanReadScrubber(Tool):

    """
    NCBI human read scrubbing tool, also called HRRT or human read removal tool.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the HRRT.
        :param camel: Camel instance
        """
        super().__init__('HRRT', '2.2.1', camel)

    def _execute_tool(self) -> None:
        """
        Runs the HRRT.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        if 'FASTQ_SINGLE_GUNZIP' not in self._tool_inputs or len(self._tool_inputs['FASTQ_SINGLE_GUNZIP']) == 0:
            raise ValueError("Required FASTQ input file is missing for human read scrubber.")

    def __build_command(self) -> None:
        """
        Builds the command line call to execute HRRT.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join(self._build_options(excluded_parameters=['interleaved'])),
            self._parameters['interleaved'].option if self._parameters['interleaved'].value == 'true' else '',
            '-i', str(self._tool_inputs['FASTQ_SINGLE_GUNZIP'][0].path)])

    def _check_command_output(self) -> None:
        """
        Checks if the command output is valid.
        :return: None
        """
        if not self._command.returncode == 0:
            raise ToolExecutionError(f"Error executing {self.name}: {self._command.stderr}")

    def __set_output(self) -> None:
        """
        Set the output of HRRT.
        :return: None
        """
        self._tool_outputs['FASTQ_SCRUBBED'] = [ToolIOFile(Path(self._parameters['outputfile'].value))]
