import os

from app.error.toolexecutionerror import ToolExecutionError
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.io.tooliofile import ToolIOFile
from app.tools.samtools.samtools import Samtools
import logging


class SamtoolsMerge(Samtools):
    """
    Merges bam/sam files into one.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        :return: None
        """
        super(SamtoolsMerge, self).__init__('samtools merge', '1.3.1', camel)

    def _check_input(self):
        """
        Checks the validity of input file types (SAM/BAM) and sets the file type accordingly.
        :return: None
        """
        if 'BAM' not in self._tool_inputs and 'SAM' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No BAM or SAM file given as input.")
        elif 'BAM' in self._tool_inputs and 'SAM' in self._tool_inputs:
            raise InvalidInputSpecificationError("BAM and SAM files given as input; tool can only accept one type.")
        elif 'BAM' in self._tool_inputs:
            self.__input_file_type = 'BAM'
        elif 'SAM' in self._tool_inputs:
            self.__input_file_type = 'SAM'

        super(Samtools, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__output_file_type = self.__determine_output_file_type()
        self.__build_command()
        self._execute_command()
        self._check_stderr()
        self.__set_output()


    def __build_command(self):
        """
        Builds the command for this tool.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join(self._build_options(excluded_parameters = ('output_filename'))),
            self._parameters['output_filename'].value,
            ' '.join([file.path for file in self._tool_inputs[self.__input_file_type]])
            ])


    def _check_stderr(self):
        """
        Validates the stderr.
        :return: None
        TODO
        """

        super(SamtoolsMerge, self)._check_stderr()


    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """
        output_file_path = os.path.join(self._folder, self._parameters['output_filename'].value)
        self._tool_outputs[self.__output_file_type] = [ToolIOFile(output_file_path)]


    def __determine_output_file_type(self):
        """
        Determines the output file type to use based on output_filename.
        :return: filetype string ("BAM" or "SAM") 
        """
        if self._parameters['output_filename'].value.split(".")[-1].lower() == "bam":
            filetype = "BAM"
        elif self._parameters['output_filename'].value.split(".")[-1].lower() == "sam":
            filetype = "SAM"
        else:
            filetype = "BAM"
            logging.info("Output file format not BAM or SAM. Assuming BAM and continuing.")
        return filetype

    def _check_command_output(self):
        pass