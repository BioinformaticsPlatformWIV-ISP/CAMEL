from pathlib import Path

from camel.app.camel import Camel
from camel.app.tools.picard.picard import Picard
from camel.app.io.tooliofile import ToolIOFile

class CollectRawWgsMetrics(Picard):
    """
    Class for picard CollectRawWgsMetrics function. Collect metrics about coverage and performance of WGS experiments.
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard CollectRawWgsMetrics', '2.23.3', camel)

        self._function_name = 'CollectRawWgsMetrics'

        # individual files of different types that is required: e.g, FASTA_REF
        self._required_inputs = ['FASTA_REF']

    def _set_input(self) -> None:
        """
        Set input specification
        :return: None
        """
        if 'COVERAGE_INTERVALS' in self._tool_inputs:
            self._input_string += f" INTERVALS={self._tool_inputs['COVERAGE_INTERVALS'][0].path}"

        if 'FASTA_REF' in self._tool_inputs:
            self._input_string += f" R={self._tool_inputs['FASTA_REF'][0].path}"

    def _set_output(self) -> None:
        """
        Set the output specification. Overrides method in parent class.
        :return: None
        """
        self._tool_outputs['TXT_metrics'] = [ToolIOFile(Path(self._folder) / self._parameters['output'].value)]