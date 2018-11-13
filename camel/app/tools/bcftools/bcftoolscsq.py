import os

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class BcftoolsCsq(Tool):
    """
    Bcftools csq is a Haplotype-aware consequence caller.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL
        """
        super().__init__('bcftools csq', '1.9', camel)

    def _check_input(self):
        """
        Checks if the provided input is valid.
        :return: None
        """
        if not any(key in self._tool_inputs for key in ('VCF', 'VCF_GZ')):
            raise InvalidInputSpecificationError("VCF(_GZ) input is required.")
        super()._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self._tool_outputs['VCF'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output_filename'].value))]

    @property
    def _input_key(self):
        """
        Returns the input key.
        :return: None
        """
        return 'VCF' if 'VCF' in self._tool_inputs else 'VCF_GZ'

    def __build_command(self):
        """
        Builds the command line command.
        :return: None
        """
        parts = [
            self._tool_command,
            ' '.join(self._build_options()),
            self._tool_inputs[self._input_key][0].path]
        if 'GFF' in self._tool_inputs:
            parts.insert(3, '--gff-annot {}'.format(self._tool_inputs['GFF'][0].path))
        if 'FASTA' in self._tool_inputs:
            parts.insert(3, '--fasta-ref {}'.format(self._tool_inputs['FASTA'][0].path))
        self._command.command = ' '.join(parts)

    def _check_command_output(self):
        """
        Checks if the command executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Error executing {}:\n{}".format(self._name, self._command.stderr))
