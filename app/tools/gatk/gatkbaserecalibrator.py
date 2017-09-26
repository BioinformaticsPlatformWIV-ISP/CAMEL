import os

from app.io.tooliofile import ToolIOFile
from app.tools.gatk.gatk import GATK
from app.io.tooliodb import ToolIODb

import logging


class GATKBaseRecalibrator(GATK):
    """
    Class for GATK BaseRecalibrator tool.
    """

    def __init__(self, camel):
        """
        Initialize GATKBaseRecalibrator tool.
        :param camel: Camel instance
        :return: None
        """
        super(GATKBaseRecalibrator, self).__init__('gatk BaseRecalibrator', '3.7', camel)

        self._function_name = 'BaseRecalibrator'
        self._required_inputs = ['BAM','FASTA_REF']
        self._output_type = 'TXT_RecalibrationTable'
        self.__snps_known_sites_path = ''
        self.__indels_known_sites_path = ''

    def _set_input(self):
        """
        Set the input specification in the input_string
        :return: None
        """

        # set input BAM
        self._input_string += "-I {} ".format(self._tool_inputs['BAM'][0].path)

        # set reference genome, known snps and known indels
        self._input_string += "-R {} ".format(self._tool_inputs['FASTA_REF'][0].path)


        if 'VCF_KNOWN_SNPS' in self._tool_inputs:
            self._input_string += "-knownSites {} ".format(self._tool_inputs['VCF_KNOWN_SNPS'][0].path)

        if 'VCF_KNOWN_INDELS' in self._tool_inputs:
            self._input_string += "-knownSites {} ".format(self._tool_inputs['VCF_KNOWN_INDELS'][0].path)

        if 'TXT_intervals' in self._tool_inputs:
            self._input_string += "-L {} ".format(self._tool_inputs['TXT_intervals'][0].path)

    def _set_output(self):
        """
        Set the output specification in the output_string
        :return: None
        """
        self._tool_outputs['TXT_RecalibrationTable'] = [
            ToolIOFile(os.path.join(self._folder, self._parameters['recal_table_output'].value))]
