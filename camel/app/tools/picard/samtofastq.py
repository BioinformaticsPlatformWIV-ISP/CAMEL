import os

from camel.app.camel import Camel
from camel.app.tools.picard.picard import Picard
from camel.app.io.tooliofile import ToolIOFile


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

    def _set_output(self) -> None:
        """
        Set the output specification
        :return: None
        """
        self._tool_outputs['FASTQ'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output'].value))]