import os

from app.io.tooliofile import ToolIOFile
from app.tools.gatk.gatk import GATK
from app.io.tooliodb import ToolIODb
from app.error.invalidparametererror import InvalidParameterError
import logging


class GATKAnalyzeCovariates(GATK):
    """
    Class for the GATK AnalyzeCovariates tool.
    """
    def __init__(self,camel):
        """
        Initialize GATKAnalyzeCovariates tool.
        :param camel: Camel instance
        :return: None
        """

        super(GATKAnalyzeCovariates, self).__init__('gatk AnalyzeCovariates', '3.7', camel)

        self._function_name = 'AnalyzeCovariates'
        self._required_inputs = ['TABLE_BEFORE', 'TABLE_AFTER']


    def _set_input(self):
        """
        Set the input specification in the input_string: 
        - before and after covariates tables (and optional BQSR table)
        - reference fasta (default or superseded).
        :return: None
        """

        # set before and after covariates tables
        if 'TABLE_BEFORE' in self._tool_inputs:
            self._input_string += "-before {} ".format(self._tool_inputs['TABLE_BEFORE'][0].path)
        if 'TABLE_AFTER' in self._tool_inputs:
            self._input_string += "-after {} ".format(self._tool_inputs['TABLE_AFTER'][0].path)
        if 'BQSR' in self._tool_inputs:
            self._input_string += "-BQSR {} ".format(self._tool_inputs['BQSR'][0].path)

        # set reference genome
        if 'FASTA_REF' in self._tool_inputs:
            self._input_string += "-R {} ".format(self._tool_inputs['FASTA_REF'][0].path)
        else:
            # set default
            self.__fasta_ref = ToolIODb('broad_b37_human_Genome_1K_v37')
            self._input_string += "-R {} ".format(self.__fasta_ref)
            logging.info("Setting fasta reference to default: {}".format(self.__fasta_ref))


    def _set_output(self):
        """
        Set the output specifications in the ouptut_string: 
        - plots pdf file, 
        - data csv file (optional)
        Supersedes _set_output in GATK class.
        :return: None
        """

        self._tool_outputs['PDF'] = [ToolIOFile(os.path.join(self._folder, self._parameters['pdf_output'].value))]

        if self._parameters['write_csv_output'].value == 'True':
            if 'csv_output' in self._parameters:
                self._tool_outputs['CSV'] = [ToolIOFile(os.path.join(self._folder, self._parameters['csv_output'].value))]



    def _set_specific_parameters(self):
        """
        Set excluded_parameters for build_option.
        - write_csv_output (internal)
        - csv_output if not required by input.
        Overrides the set_specific_parameters fct of GATK class.
        :return: None
        """
        super(GATKAnalyzeCovariates, self)._set_specific_parameters()
        self._specific_parameters.append('write_csv_output')

        if self._parameters['write_csv_output'].value == 'False':
            self._specific_parameters.append('csv_output')


    def _check_parameters(self):
        """
        Check that parameters make sense and that mandatory parameters were specified.
        :return: None
        """

        super(GATKAnalyzeCovariates, self)._check_parameters()

        if self._parameters['write_csv_output'].value not in ('False','True'):
            raise InvalidParameterError("Unrecognized boolean value: {}. Value should be 'True' or 'False'".format(
                self._parameters['write_csv_output'].value))