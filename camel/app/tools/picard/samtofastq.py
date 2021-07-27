import os

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.picard.picard import Picard


class SamToFastq(Picard):

    """
    Class for Picard FastqToSam function
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard SamToFastq', '2.23.3', camel)
        self._function_name = 'SamToFastq'
        self._main_input = []

    def _check_input(self) -> None:
        """
        Check the main inputs and set _main_input
        _main_input can be either SAM or BAM, not both
        :return: None
        """
        super(Picard, self)._check_input()

        for input_format in self._main_inputs:
            if input_format in self._tool_inputs:
                self._main_input.append(input_format)

        if len(self._main_input) != 1:
            raise InvalidInputSpecificationError("Picard MarkDuplicates requires one SAM or BAM file")

    def _set_input(self) -> None:
        """
        Set the input specification
        :return: None
        """
        # input format can be either SAM or BAM
        main_input = self._main_input[0]

        self._input_string = f"I={self._tool_inputs[main_input][0].path}"

    def _set_output(self) -> None:
        """
        Set the output specification
        :return: None
        """
        self._tool_outputs['FASTQ'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output'].value))]