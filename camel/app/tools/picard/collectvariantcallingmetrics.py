from camel.app.camel import Camel
from camel.app.tools.picard.picard import Picard

class CollectVariantCallingMetrics(Picard):
    """
    Class for picard CollectVariantCallingMetrics function
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard CollectVariantCallingMetrics', '2.23.3', camel)

        self._function_name = 'CollectVariantCallingMetrics'
        self._supported_inputs = ['VCF', 'VCF_dbsnp', 'DICT_GENOME', 'EVALUATION_INTERVALS']
        self._required_inputs = ['VCF', 'VCF_dbsnp']

    def _check_input(self) -> None:
        """
        Check and set the input specification.  Overrides method in parent class.
        :return: None
        """
        super(Picard, self)._check_input()

        self._set_input()

    def _set_input(self) -> None:
        """
        Set the input specification. Overrides method in parent class.
        :return: None
        """
        self._input_string += f" INPUT={self._tool_inputs['VCF'][0].path}"

        self._input_string += f" DBSNP={self._tool_inputs['VCF_dbsnp'][0].path}"

        if 'DICT_GENOME' in self._tool_inputs:
            self._input_string += f" SEQUENCE_DICTIONARY={self._tool_inputs['DICT_GENOME'][0].path}"

        if 'EVALUATION_INTERVALS' in self._tool_inputs:
            self._input_string += f" TARGET_INTERVALS={self._tool_inputs['EVALUATION_INTERVALS'][0].path}"

    def _set_output(self) -> None:
            """
            Set the output specification. Overrides method in parent class.
            :return: None
            """
            self._tool_outputs['TXT_report'] = [
                ToolIOFile(os.path.join(self._folder, self._parameters['output'].value, ".variant_calling_summary_metrics")),
                ToolIOFile(os.path.join(self._folder, self._parameters['output'].value, ".variant_calling_detail_metrics"))
            ]