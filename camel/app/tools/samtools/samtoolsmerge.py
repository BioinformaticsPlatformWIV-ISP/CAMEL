import os

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.samtools.samtools import Samtools
import logging


class SamtoolsMerge(Samtools):
    """
    ==============
    SamtoolsMerge 1.9.
    ==============
    Merges bam/sam files into one.

    required inputs:
    ----------------
    "BAM" / "SAM":  At least one bam or sam file. All files should be of same type.

    Output:
    -------
    "BAM" / "SAM":  one bam or sam file. Output filetype depends on output file extension; bam assumed if no valid
                    extension is found.

    Mandatory parameters:
    ---------------------
    - output_filename
                    default value:  merged.bam
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        :return: None
        """
        super(SamtoolsMerge, self).__init__('samtools merge', '1.9', camel)

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
            ' '.join(self._build_options(excluded_parameters={'output_filename'})),
            self._parameters['output_filename'].value,
            ' '.join([file.path for file in self._tool_inputs[self.__input_file_type]])
        ])

    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """
        output_file_type = self.__determine_output_file_type()
        output_file_path = os.path.join(self._folder, self._parameters['output_filename'].value)
        self._tool_outputs[output_file_type] = [ToolIOFile(output_file_path)]

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
        """
        Checks the command output.
        Supersedes function in Tool class because warnings printed to stderr can cause false abort.
        """
        self._check_stderr()

        if self._command.returncode != 0:
            raise ToolExecutionError("Command execution failed (Exit code: {})".format(self._command.returncode))
