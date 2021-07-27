import os

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.picard.picard import Picard


class ValidateSamFile(Picard):
    """
    Class for Picard ValidateSamFile function
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard ValidateSamFile', '2.23.3', camel)
        self._extra_inputs = ["FASTA_REF"]

    def _set_input(self) -> None:
        """
        Set input specification
        """
        super(ValidateSamFile, self)._set_input()

        if 'FASTA_REF' in self._tool_inputs:
            self._input_string += f'R={self._tool_inputs["FASTA_REF"][0].path} '

    def _set_output(self) -> None:
        """
        Set the output specification, this function handles a TXT_report as output
        :return: None
        """
        self._tool_outputs['TXT_report'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output'].value))]