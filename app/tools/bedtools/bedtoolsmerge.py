import os

from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.tools.bedtools.bedtools import Bedtools
from app.io.tooliofile import ToolIOFile


class BedtoolsMerge(Bedtools):
    """
    Tool class for Bedtools merge function.
    """

    def __init__(self, camel, tool_name='bedtools merge', version='2.25.0'):
        """
        Initialize a samtools tool.
        :param tool_name: Tool name
        :param version: Tool version
        :param camel: Camel instance
        :return: None
        """
        super(BedtoolsMerge, self).__init__(tool_name, version, camel)
        # self._required_inputs = ['BED']
        self.__input_type = ""

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self._check_input()
        if self.__input_type == "BAM":
            self.update_parameters(use_bed_output=True)

        self.__set_output()
        self.__build_command()
        self._execute_command()

    def __build_command(self):
        """
        Builds the command with input, options and output strings.
        :return: None
        """

        input_string = "-i {} ".format(self._tool_inputs[self.__input_type][0].path)

        build_options = ' '.join(self._build_options(excluded_parameters='output_filename'))

        output_string = '> ' + self._parameters['output_filename'].value

        self._command.command = ' '.join([
            self._tool_command,
            build_options,
            input_string,
            output_string])

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        self._check_required_inputs()
        super(BedtoolsMerge, self)._check_input()

    def _check_required_inputs(self):
        """
        Checks that required input is present and sets the input type.
        Supersedes the same function in the Bedtools class.
        Inputs must be one of either BAM or BED.
        :return: None
        """

        if len(self._tool_inputs) != 1:
            raise InvalidInputSpecificationError(
                "{} input file(s) specified. Bedtools merge takes exactly ONE input (BAM or BED).".format(
                    len(self._tool_inputs)))
        elif "BAM" in self._tool_inputs:
            self.__input_type = "BAM"
        elif "BED" in self._tool_inputs:
            self.__input_type = "BED"
        else:
            raise InvalidInputSpecificationError(
                "Input file specified with wrong file type ({}). Accepted types are BAM or BED.".format(
                    self._tool_inputs.keys()[0]))

    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """

        self._tool_outputs['BED'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output_filename'].value))]
