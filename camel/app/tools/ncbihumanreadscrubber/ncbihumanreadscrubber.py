import re

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
        self._parse_stderr()
        self.__set_output()

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        if 'FASTQ_SINGLE_GUNZIP' not in self._tool_inputs or len(self._tool_inputs['FASTQ_SINGLE_GUNZIP']) == 0:
            raise ValueError("Required FASTQ input file is missing for human read scrubber.")
        super()._check_input()

    def __build_command(self) -> None:
        """
        Builds the command line call to execute HRRT.
        Export_human_reads and outputfile_removed linked, adds args -r -u and path if export_human_reads = true
        :return: None
        """
        dir_temp = self._camel.config['temp_dir']
        self._command.command = ' '.join([
            f'export TMPDIR={dir_temp};',
            self._tool_command,
            *self._build_options(excluded_parameters=['interleaved', 'export_human_reads', 'outputfile_removed']),
            self._parameters['interleaved'].option if self._parameters['interleaved'].value == 'true' else '',
            self._parameters['export_human_reads'].option if self._parameters['export_human_reads'].value == 'true' else '',
            self._parameters['outputfile_removed'].option + ' ' + self._parameters['outputfile_removed'].value if self._parameters['export_human_reads'].value == 'true' else '',
            '-i', str(self._tool_inputs['FASTQ_SINGLE_GUNZIP'][0].path)
            ])

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
        path_out = self.folder / self._parameters['outputfile'].value
        path_removed = self.folder / self._parameters['outputfile_removed'].value if self._parameters['export_human_reads'].value == 'true' else []
        self._tool_outputs['FASTQ_SCRUBBED'] = [ToolIOFile(path_out)]
        if self._informs.get('statistics').get('count_removed') != 0 and self._parameters['export_human_reads'].value == 'true':
            self._tool_outputs['FASTQ_REMOVED'] = [ToolIOFile(path_removed)]
        else:
            self.tool_outputs['FASTQ_REMOVED'] = []

    def _parse_stderr(self) -> None:
        """
        Parses the command's stderr to retrieve the statistics about how many reads/contigs were removed.
        :return: None
        """
        count_removed = None
        count_total = None
        for line in self._command.stderr.split('\n'):
            # Define the regular expression pattern
            pattern_reads_removed = r'^(\d+)\s+spot\(s\) masked or removed\.$'
            pattern_reads_total = r'total spot count: (\d+)'

            # Try to match the pattern in the current line
            if count_removed is None and re.match(pattern_reads_removed, line):
                # Extract the matched integer
                count_removed = int((re.match(pattern_reads_removed, line)).group(1))
            elif count_total is None and re.search(pattern_reads_total, line):
                # Extract the matched integer
                count_total = int((re.search(pattern_reads_total, line)).group(1))
            elif count_removed is not None and count_total is not None:
                break
        self._informs['statistics'] = {'count_removed': count_removed, 'count_total': count_total}
