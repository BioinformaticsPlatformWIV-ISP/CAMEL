import os

from app.io.tooliofile import ToolIOFile
from app.tools.gatk.gatk import GATK


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
    'FASTA_REF':        FASTA file containing the reference genome. If not specified, db default is used.
    
    Optional input:
    ---------------
    'BQSR':             Base quality score recalibration table on extra data. GATK text format.
    
    Output:
    -------
    'PDF':              pdf document with graphs of recalibration.
    'CSV':              csv with recalibration covariates data used for generating the pdf output. 
    
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
        self._required_inputs = ['TXT_TABLE_BEFORE', 'TXT_TABLE_AFTER', 'FASTA_REF']

    def _set_input(self):
        """
        Set the input specification in the input_string: 
        - before and after covariates tables (and optional BQSR table)
        - reference fasta (default or superseded).
        :return: None
        """
        self._input_string += "-before {} ".format(self._tool_inputs['TXT_TABLE_BEFORE'][0].path)
        self._input_string += "-after {} ".format(self._tool_inputs['TXT_TABLE_AFTER'][0].path)

        self._input_string += "-R {} ".format(self._tool_inputs['FASTA_REF'][0].path)

        if 'BQSR' in self._tool_inputs:
            self._input_string += "-BQSR {} ".format(self._tool_inputs['BQSR'][0].path)

    def _set_output(self):
        """
        Set the output specifications in the ouptut_string: 
        - plots pdf file, 
        - data csv file (optional)
        Supersedes _set_output in GATK class.
        :return: None
        """
        self._tool_outputs['PDF'] = [ToolIOFile(os.path.join(self._folder, self._parameters['pdf_output'].value))]

        self._tool_outputs['CSV'] = [ToolIOFile(os.path.join(self._folder, self._parameters['csv_output'].value))]

    def _check_parameters(self):
        """
        Check that parameters make sense and that mandatory parameters were specified.
        :return: None
        """

        super(GATKAnalyzeCovariates, self)._check_parameters()
