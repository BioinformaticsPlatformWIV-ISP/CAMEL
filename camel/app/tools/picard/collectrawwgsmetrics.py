from camel.app.camel import Camel
from camel.app.tools.picard.picard import Picard


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
        self._required_inputs = ['BAM', 'SAM', 'FASTA_REF']
        self._output_type = 'TXT_metrics'

    def _set_input(self) -> None:
        """
        Set input specification
        :return: None
        """
        super(CollectRawWgsMetrics, self)._set_input()

        if 'COVERAGE_INTERVALS' in self._tool_inputs:
            self._input_string += f" INTERVALS={self._tool_inputs['COVERAGE_INTERVALS'][0].path}"