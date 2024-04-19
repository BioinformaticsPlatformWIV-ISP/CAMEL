import re

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.invalidparametererror import InvalidParameterError
from camel.app.loggers import logger
from camel.app.tools.gatk.gatk import GATK


class GATKSelectVariants(GATK):

    """
    Class for GATK SelectVariants function
    """

    def __init__(self, camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('gatk SelectVariants', '3.7', camel)
        self._specific_parameters = ['selectTypeToExlcude', 'selectTypeToInclude', 'select']
        self._required_inputs = ['VCF']
        self._output_type = 'VCF'

    def _check_parameters(self):
        """
        Check the parameters set to run this tool
        :return: None
        """
        if 'select' in self._parameters:
            select_exp = self._parameters['select'].value
            if not select_exp:
                raise InvalidParameterError(
                    "Selection criterion of 'select' option is required.")
            if re.search('\s', select_exp):
                raise InvalidParameterError(
                    "No space in the 'select' option JEXL expression. e.g., 'DQ>100', not 'DQ > 100'. Expression specified {}.".format(select_exp))
        else:
            logger.warning("GATK SelectVariants run without specify variant selection criterion.")

    def _check_input(self):
        """
        Check input for a tool and prepare command line parameters for input
        :return: None
        """
        super(GATKSelectVariants, self)._check_input()

        if 'VCF_concordance' in self._tool_inputs and len(self._tool_inputs['VCF_concordance']) > 1:
            raise InvalidInputSpecificationError('SelectVariant support only ONE concordance VCF file.')

        if 'VCF_discordance' in self._tool_inputs and len(self._tool_inputs['VCF_discordance']) > 1:
            raise InvalidInputSpecificationError('SelectVariant support only ONE discordance VCF file.')

    def _set_input(self):
        """
        Set the input specification
        :return: None
        """
        super(GATKSelectVariants, self)._set_input()

        # To extract only certain samples
        if 'SAMPLES' in self._tool_inputs:
            for sample_name in self._tool_inputs['SAMPLES']:
                self._input_string += '--sample_name {} '.format(sample_name.value)

        # To exclude certain samples
        if 'XL_SAMPLES' in self._tool_inputs:
            for sample_name in self._tool_inputs['XL_SAMPLES']:
                self._input_string += '--exclude_sample_name {} '.format(sample_name.value)

        # Variant comparison (with those reported in another VCF file)
        # - concordance only
        if 'VCF_concordance' in self._tool_inputs:
            self._input_string += '--concordance {} '.format(self._tool_inputs['VCF_concordance'][0].path)
        # - discordance only
        if 'VCF_discordance' in self._tool_inputs:
            self._input_string += '--discordance {} '.format(self._tool_inputs['VCF_discordance'][0].path)

    def _set_specific_parameters(self):
        """
        Set specific parameters that need special handling
        :return: None
        """
        if 'selectTypeToExclude' in self._parameters:
            option = self._parameters['selectTypeToExclude'].option
            value = self._parameters['selectTypeToExclude'].value
            excluded_types = value.split(",")
            for ex_type in excluded_types:
                self._option_string += "{} {!r} ".format(option, ex_type)

        if 'selectTypeToInclude' in self._parameters:
            option = self._parameters['selectTypeToInclude'].option
            value = self._parameters['selectTypeToInclude'].value
            selected_types = value.split(",")
            for se_type in selected_types:
                self._option_string += "{} {!r} ".format(option, se_type)

        if 'select' in self._parameters:
            select_opt = self._parameters['select']
            self._option_string += "{} {!r} ".format(select_opt.option, select_opt.value)
