import os

from app.io.tooliofile import ToolIOFile
from app.tools.gatk.gatk import GATK
from app.io.tooliodb import ToolIODb
from app.error.invalidparametererror import InvalidParameterError
import logging


class GATKAnalyzeCovariates(GATK):
    """
    ==========================
    GATK AnalyzeCovariates 3.7
    ==========================
    Analyses and generates plots about covariates analysed during Base Quality Score Recalibration step and used to recalibrate bases.
    web: https://software.broadinstitute.org/gatk/gatkdocs/3.5-0/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php
    
    Required inputs:
    ----------------
    'TXT_TABLE_BEFORE': Base quality score recalibration table on original data. GATK text format.
    'TXT_TABLE_AFTER':  Base quality score recalibration table on recalibrated data. GATK text format.
    
    Optional input:
    ---------------
    'BQSR':             Base quality score recalibration table on extra data. GATK text format.
    'FASTA_REF':        FASTA file containing the reference genome. If not specified, db default is used.
     
    Output:
    -------
    'PDF':              pdf document with graphs of recalibration.
    
    Optional output:
    ---------------
    'CSV':              csv with recalibration covariates data used for generating the pdf output. 
                        Generated if 'csv_output' parameter set to 'True'
    
    Mandatory parameters:
    ---------------------
    - pdf_output       Default value: 'recal_QC_plots.pdf'
    """

    def __init__(self, camel):
        """
        Initialize GATKAnalyzeCovariates tool.
        :param camel: Camel instance
        :return: None
        """

        super(GATKAnalyzeCovariates, self).__init__('gatk AnalyzeCovariates', '3.7', camel)

        self._function_name = 'AnalyzeCovariates'
        self._required_inputs = ['TXT_TABLE_BEFORE', 'TXT_TABLE_AFTER']

    def _set_input(self):
        """
        Set the input specification in the input_string: 
        - before and after covariates tables (and optional BQSR table)
        - reference fasta (default or superseded).
        :return: None
        """

        # set before and after covariates tables
        if 'TXT_TABLE_BEFORE' in self._tool_inputs:
            self._input_string += "-before {} ".format(self._tool_inputs['TXT_TABLE_BEFORE'][0].path)
        if 'TXT_TABLE_AFTER' in self._tool_inputs:
            self._input_string += "-after {} ".format(self._tool_inputs['TXT_TABLE_AFTER'][0].path)
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
                self._tool_outputs['CSV'] = [
                    ToolIOFile(os.path.join(self._folder, self._parameters['csv_output'].value))]

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

        if self._parameters['write_csv_output'].value not in ('False', 'True'):
            raise InvalidParameterError("Unrecognized boolean value: {}. Value should be 'True' or 'False'".format(
                self._parameters['write_csv_output'].value))
