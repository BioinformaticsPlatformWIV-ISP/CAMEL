from pathlib import Path

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
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
        self._main_inputs = ['VCF']
        self._extra_inputs = ['VCF_dbsnp', 'DICT_GENOME', 'EVALUATION_INTERVALS']

    def _check_input(self) -> None:
        """
        Check and set the input specification.  Overrides method in parent class.
        :return: None
        """
        super(Picard, self)._check_input()

        if 'VCF_dbsnp' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Picard CollectVariantCallingMetrics: input file VCF_dbSNP is not specified")

    def _set_input(self) -> None:
        """
        Set the input specification. Overrides method in parent class.
        :return: None
        """
        super(CollectVariantCallingMetrics, self)._set_input()

        self._input_string += f"DBSNP={self._tool_inputs['VCF_dbsnp'][0].path} "

        if 'DICT_GENOME' in self._tool_inputs:
            self._input_string += f"SEQUENCE_DICTIONARY={self._tool_inputs['DICT_GENOME'][0].path} "

        if 'EVALUATION_INTERVALS' in self._tool_inputs:
            self._input_string += f"TARGET_INTERVALS={self._tool_inputs['EVALUATION_INTERVALS'][0].path} "

    def _set_output(self) -> None:
            """
            Set the output specification. Overrides method in parent class.
            :return: None
            """
            self._tool_outputs['TXT_report'] = [
                ToolIOFile(Path(self._folder) / f"{self._parameters['output_prefix'].value}.variant_calling_summary_metrics"),
                ToolIOFile(Path(self._folder) / f"{self._parameters['output_prefix'].value}.variant_calling_detail_metrics")
            ]