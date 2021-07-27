from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.picard.picard import Picard


class CollectWgsMetrics(Picard):
    """
    Class for picard CollectWgsMetrics function. Collect metrics about coverage and performance of WGS experiments.
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard CollectWgsMetrics', '2.23.3', camel)

        self._function_name = 'CollectWgsMetrics'
        self._extra_inputs = ['FASTA_REF']

    def _check_input(self) -> None:
        """
        Check input specification
        :return: None
        """
        super(CollectWgsMetrics, self)._check_input()

        if 'FASTA_REF' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Picard CollectWgsMetrics: input file FASTA_REF is not defined")

    def _set_input(self) -> None:
        """
        Set input specification
        :return: None
        """
        super(CollectWgsMetrics, self)._set_input()

        self._input_string += f" R={self._tool_inputs['FASTA_REF'][0].path}"

        if 'COVERAGE_INTERVALS' in self._tool_inputs:
            self._input_string += f" INTERVALS={self._tool_inputs['COVERAGE_INTERVALS'][0].path}"

    def _set_output(self) -> None:
        """
        Set the output specification
        :return: None
        """
        self._tool_outputs['TXT_metrics'] = [ToolIOFile(Path(self._folder) / self._parameters['output'].value)]