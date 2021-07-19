from pathlib import Path

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.picard.picard import Picard


class CalculateReadGroupChecksum(Picard):

    """
    ==============================
    Picard CalculateReadGroupChecksum 2.23.3
    ==============================
    Creates a hash code based on the read groups (RG)

    Output:
    -------
    'TXT_checksum'  ToolIOFile object. File to which the hash code will be written
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard CalculateReadGroupChecksum', '2.23.3', camel)

        self._function_name = 'CalculateReadGroupChecksum'

    def _set_output(self) -> None:
        """
        Set the output specification
        Overrides method from the parent class
        :return: None
        """
        self._tool_outputs['TXT_checksum'] = [ToolIOFile(Path(self._folder) / self._parameters['output'].value)]