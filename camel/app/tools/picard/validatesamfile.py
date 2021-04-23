import os

from camel.app.camel import Camel
from camel.app.tools.picard.picard import Picard
from camel.app.io.tooliofile import ToolIOFile


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

    def _set_output(self) -> None:
        """
        Set the output specification, this default function handles only one BAM file as output
        :return: None
        """
        self._tool_outputs['TXT_report'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output'].value))]