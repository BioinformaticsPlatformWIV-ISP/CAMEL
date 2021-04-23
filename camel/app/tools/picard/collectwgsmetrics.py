from camel.app.camel import Camel
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

        self._supported_inputs = ['SAM', 'BAM']
        # individual files of different types that is required: e.g, FASTA_REF
        self._required_inputs = ['FASTA_REF']

    def _check_input(self) -> None:
        """
        Check and set the input specification
        :return: None
        """
        super(Picard, self)._check_input()

        self._set_input()

    def _set_input(self) -> None:
        """
        Set input specification
        :return: None
        """

        self._input_string = f"--INPUT {self._tool_inputs['BAM'][0].path} "

        self._input_string = f"--REFERENCE_SEQUENCE {self._tool_inputs['FASTA_REF'][0].path} "

        if 'COVERAGE_INTERVALS' in self._tool_inputs:
            self._input_string = f"--INTERVALS {self._tool_inputs['COVERAGE_INTERVALS'][0].path} "

    def _set_output(self) -> None:
        """
        Set the output specification
        :return: None
        """
        self._tool_outputs['TXT_metrics'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output'].value))]