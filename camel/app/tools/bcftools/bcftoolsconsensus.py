import os

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class BcftoolsConsensus(Tool):
    """
    Create consensus sequence by applying VCF variants to a reference fasta file.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('bcftools consensus', '1.9', camel)

    def _check_input(self):
        """
        Checks if the input is valid.
        :return: None
        """
        if 'VCF_GZ' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No variants input found.")
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No reference FASTA input found.")
        super(BcftoolsConsensus, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def __build_command(self):
        """
        Builds the command for this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            '--fasta-ref {}'.format(self._tool_inputs['FASTA'][0].path),
            self._tool_inputs['VCF_GZ'][0].path,
            ' '.join(self._build_options())
        ])

    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['FASTA'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output_filename'].value))]

    def _check_command_output(self):
        """
        Checks if the command executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Tool did not run successfully")
