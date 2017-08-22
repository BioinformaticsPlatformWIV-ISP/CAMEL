import os

from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.error.invalidparametererror import InvalidParameterError
from app.error.toolexecutionerror import ToolExecutionError
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class PhenixVcfs2Fasta(Tool):
    """
    Converts the VCF files to a SNP matrix in FASTA format.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(PhenixVcfs2Fasta, self).__init__('PHEnix Vcfs2Fasta', '1.2', camel)

    def _check_input(self):
        """
        Checks if the input is valid.
        VALID inputs:
        - DIR + regex (opt)
        - List of VCF files
        :return: None
        """
        if not any (key in self._tool_inputs for key in ('DIR_VCF', 'VCF')):
            raise InvalidInputSpecificationError("Either DIR_VCF or VCF is required as input.")
        super(PhenixVcfs2Fasta, self)._check_input()

    def _check_parameters(self):
        """
        Checks if the given parameters are valid.
        :return: None
        """
        if 'regexp' in self._parameters and 'DIR_VCF' not in self._tool_inputs:
            raise InvalidParameterError("Regexp parameter is only used when the input is a directory")
        super(PhenixVcfs2Fasta, self)._check_parameters()

    def _execute_tool(self):
        """
        Runs this tool.
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()
        self.__set_informs()

    def __build_command(self):
        """
        Builds the command.
        :return: None
        """
        command_parts = ['. $VIRTUALENV;', self._tool_command]

        if 'DIR_VCF' in self._tool_inputs:
            command_parts.append('--directory {}'.format(self._tool_inputs['DIR_VCF'][0].path))
        else:
            command_parts.append('--input {}'.format(' '.join([f.path for f in self._tool_inputs['VCF']])))
            command_parts.append(' '.join(self._build_options()))
        self._command.command = ' '.join(command_parts)

    def __set_output(self):
        """
        Sets the tool output.
        :return: None
        """
        self._tool_outputs['FASTA'] = [ToolIOFile(os.path.join(
            self._folder, self._parameters['output_filename'].value))]

    def __set_informs(self):
        """
        Sets the tool informs.
        :return: None
        """
        try:
            with open(self._tool_outputs['FASTA'][0].path) as handle:
                lines = handle.readlines()
                if len(lines) == 0:
                    raise ToolExecutionError("Output file is empty")
                self._informs['size'] = len(lines[1].strip())
        except IOError:
            raise ToolExecutionError("No output file generated.")

    def _check_command_output(self):
        """
        Checks the command output to check if the program ran correctly.
        :return: None
        """
        if 'error' in self.stderr.lower():
            raise ToolExecutionError("Error executing vcf2fasta: {}".format(self.stderr))
        if 'No VCFs found' in self.stderr:
            raise ToolExecutionError("No VCF input files found")
