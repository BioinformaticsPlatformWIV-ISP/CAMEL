import os

from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class PhenixPrepareReference(Tool):
    """
    Prepares the reference for the PHEnix SNP pipeline.
    """

    def __init__(self, camel):
        """
        Initializes the tool.
        :param camel: CAMEL instance
        """
        super(PhenixPrepareReference, self).__init__('PHEnix Prepare Reference', '1.2', camel)

    def _check_input(self):
        """
        Checks the tool input.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No FASTA input found")
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidInputSpecificationError("Only one input FASTA file supported")
        super(PhenixPrepareReference, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        if 'reference_name' in self._parameters:
            symlink_ref = os.path.join(
                self._folder, '{}.fasta'.format(FileSystemHelper.make_valid(self._parameters['reference_name'].value)))
        else:
            symlink_ref = os.path.join(self._folder, self._tool_inputs['FASTA'][0].basename)
        if os.path.islink(symlink_ref):
            os.remove(symlink_ref)
        os.symlink(self._tool_inputs['FASTA'][0].path, symlink_ref)
        self._command.command = '. $VIRTUALENV; {} --reference {} {}'.format(
            self._tool_command, symlink_ref, ' '.join(self._build_options(excluded_parameters=['reference_name'])))
        self._execute_command()
        self._tool_outputs['FASTA'] = [ToolIOFile(symlink_ref)]

    def _check_command_output(self):
        """
        Checks if the command ran correctly by analyzing the output.
        :return: None
        """
        pass
