import os

from camel.app.camel import Camel
from camel.app.tools.picard.picard import Picard
from camel.app.io.tooliofile import ToolIOFile

class MergeVCFs(Picard):
    """
    Class for picard MergeVCFs function
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard MergeVCFs', '2.23.3', camel)

        self._function_name = 'MergeVCFs'

        self._required_inputs = ['VCF']
        self._supported_inputs = ['VCF']

    def _set_input(self) -> None:
        """
        Set the input specification in the input_string
        Overrides method in parent class.
        :return: None
        """

        # set input reports
        self._input_string = ' I='
        self._input_string += ' I='.join(f.path for f in self._tool_inputs['VCF'])

    def _check_input(self) -> None:
        """
        :return: None
        """
        super(Picard, self)._check_input()

        self._set_input()

    def _set_output(self) -> None:
        """
        Set the output specification, this default function handles one VCF file as output
        :return: None
        """
        self._tool_outputs['VCF'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output'].value))]
