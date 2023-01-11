from pathlib import Path
from camel.app.camel import Camel
from camel.app.tools.tool import Tool
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile


class FastANI(Tool):
    """
    In silico taxonomic classification of Bacillus cereus group isolates using assembled genomes.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('FastANI', '1.33', camel)

    def _check_input(self) -> None:
        """
        Checks whether the required input files are specified.
        :return: None
        """
        authorized_keys = ['FASTA_Q', 'FASTA_R', 'FASTQ_R', 'FASTQ_Q',
                           'TSV_FASTA_Q', 'TSV_FASTA_R', 'TSV_FASTQ_R', 'TSV_FASTQ_Q']

        if len([key for key in self._tool_inputs if '_Q' in key]) > 1:
            raise InvalidInputSpecificationError('Please input at most one query (file or sequence)')
        if len([key for key in self._tool_inputs if '_R' in key]) > 1:
            raise InvalidInputSpecificationError('Please input at most one reference (file or sequence)')

        if len([key for key in self._tool_inputs if key in authorized_keys]) != 2:
            raise InvalidInputSpecificationError('Please check your input files')

        super()._check_input()

    def _build_command(self) -> None:
        """
        Build the command to run tool
        :return: None
        """
        fetching_query = [key for key in self._tool_inputs if '_Q' in key][0]
        query_condition = [True if 'TSV' in fetching_query else False][0]
        fetching_reference = [key for key in self._tool_inputs if '_R' in key][0]
        reference_condition = [True if 'TSV' in fetching_reference else False][0]
        if not (fetching_query and fetching_reference):
            raise InvalidInputSpecificationError(
                f'Incorrect input found: Query={fetching_query}, Ref={fetching_reference}')

        self._input_str_query = '--{} {}'.format('queryList' if query_condition else 'query',
                                                 self._tool_inputs[fetching_query][0].path)
        self._input_str_reference = '--{} {}'.format('refList' if reference_condition else 'ref',
                                                     self._tool_inputs[fetching_reference][0].path)

        self._command.command = ' '.join([self._tool_command, self._input_str_reference,
                                          self._input_str_query, *self._build_options()])

    def _check_command_output(self) -> None:
        """
        Checks command output
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Command execution failed (Exit code: {self._command.returncode})")

    def _set_output(self) -> None:
        """
        set the output file to check
        """
        output_filename = self._parameters['output_file'].value
        self._tool_outputs['TSV'] = [ToolIOFile(self.folder / Path(output_filename))]

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """

        # Building the command
        self._build_command()
        self._execute_command()

        # Collect output
        self._set_output()
