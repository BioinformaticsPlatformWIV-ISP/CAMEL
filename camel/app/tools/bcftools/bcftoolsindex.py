import os

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.samtools.samtools import Samtools


class BcftoolsIndex(Samtools):
    """
    Indexes bgzip compressed VCF files and BCF files.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super(BcftoolsIndex, self).__init__('bcftools index', '1.6', camel)
        self._input_key = None

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        if not any(key in self._tool_inputs for key in ('BCF', 'VCF_GZ')):
            raise InvalidInputSpecificationError("No input file found (BCF / VCF_GZ supported)")
        if len(self._tool_inputs) != 1:
            raise InvalidInputSpecificationError("Only one type of input is supported (VCF_GZ or BCF)")
        super(BcftoolsIndex, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self._input_key = list(self._tool_inputs.keys())[0]
        input_file_path = self.__symlink_input()
        self.__build_command(input_file_path)
        self._execute_command()
        self._check_stderr()
        self.__set_output(input_file_path)

    def __symlink_input(self):
        """
        Create a symlink for the input. This avoids cluttering the directory of the input file. This can also avoid
        errors when there are no writing permissions on the directory of the input file.
        :return: Path to symlink input
        """
        new_path = os.path.join(self._folder, self._tool_inputs[self._input_key][0].basename)
        if not os.path.isfile(new_path):
            os.symlink(self._tool_inputs[self._input_key][0].path, new_path)
        return new_path

    def __build_command(self, input_file_path):
        """
        Builds the command for this tool.
        :param input_file_path: Path to the input file
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join(self._build_options()),
            input_file_path])

    def __set_output(self, input_file_path):
        """
        Sets the output of this tool.
        :param input_file_path: Path to the input file symlink
        :return: None
        """
        self._tool_outputs[self._input_key] = [ToolIOFile(input_file_path)]
