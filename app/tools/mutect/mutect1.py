import os
import re
import logging

from app.io.tooliofile import ToolIOFile
from app.io.tooliodb import ToolIODb
from app.tools.tool import Tool
from app.error.toolexecutionerror import ToolExecutionError
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.error.invalidparametererror import  InvalidParameterError


class Mutect1(Tool):
    """
    ===========
    Mutect (v1).
    ===========
    Performs variant calling for oncology-related NGS data. 
    Mutect V1 only calls snps, not indels. For indels, use Mutect2.
    
    Required inputs:
    ----------------
    "BAM_TUMOR":        BAM file with tumour data
    
    Optional input:
    ---------------
    "BAM_NORMAL":       BAM file with normal data for tumor-normal matching.
    "FASTA_REF":        FASTA file containing the reference genome. If not specified, db default is used.
    "VCF_DBSNP":        DbSNP reference vcf file location. If not specified, db defaults is used.
    "TXT_intervals":    Intervals list to restrict search by GATK. Accelerates analysis. Bed or GATK intervals list 
    
    Output:
    -------
    "TXT_CALL_STATS": GATK Call stats text based file. Parseable by scripts or in excel sheets.
    
    Optional output:
    ---------------
    "VCF": VCF file. Generated if 'output_vcf_file' parameter set to 'True'
    
    Mandatory parameters:
    ---------------------
    - output_callstats_file
                    default value:  call_stats.txt
    """

    def __init__(self, camel):
        """
        Initialize Mutect1 tool.
        :param camel: Camel instance
        :return: None
        """
        super(Mutect1, self).__init__('Mutect1', '1.1.7', camel)

        self._required_inputs = ['BAM_TUMOR']
        self._excluded_parameters=['generate_vcf_file']

    def _execute_tool(self):
        """
        Runs Mutect1
        :return: None
        """

        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_input(self):
        """
        Check that input is valid (super method) and that required parameters are present.
        :return: None
        """
        super(Mutect1, self)._check_input()

        for input_key in self._required_inputs:
            if input_key not in self._tool_inputs:
                raise InvalidInputSpecificationError(
                    'Mutect1 required {} input is missing in tool inputs!'.format(input_key))

    def _check_parameters(self):
        """
        Checks that parameters are valid.
        :return: None 
        """
        super(Mutect1,self)._check_parameters()

        self._check_allowed_values(self._parameters['generate_vcf_file'],('True','False'))


    def _check_allowed_values(self,parameter,values):
        """
        Checks that a parameter value is one of several allowed values.
        :parameter parameter:   the parameter to test
        :parameter values:      values to test against
        :type parameter:        Parameter
        :type values:           tuple
        :return:                None
        """
        if parameter.value not in values:
            raise InvalidParameterError("Unrecognized value for parameter {}: {}. Value should be {}.".format(parameter.name,
                parameter.value, ' or '.join(values)))

    def __build_command(self):
        """
        Build the command to run the tool.
        Concatenates strings into the command string.
        :return: 
        """
        input_string = self.__create_input_string()
        options_string = ' '.join(self._build_options(excluded_parameters=self._excluded_parameters))
        self._command.command = ' '.join([self._tool_command, input_string, options_string])

    def __create_input_string(self):
        """
        Add the input specification in the input_string: 
        - DBSNP VCF file
        - reference fasta (default or superseded)
        - tumour BAM file
        - normal tissue BAM file (optional)
        - Intervals list for acceleration
        :return: Input_string
        """
        input_string = ""
        # set reference genome
        if 'FASTA_REF' in self._tool_inputs:
            input_string += "-R {} ".format(self._tool_inputs['FASTA_REF'][0].path)
        else:
            # set default
            self.__fasta_ref = ToolIODb('broad_b37_human_Genome_1K_v37')
            input_string += "-R {} ".format(self.__fasta_ref)
            logging.info("Setting fasta reference to default: {}".format(self.__fasta_ref))

        # Use intervals to restrict search if supplied.
        if 'TXT_intervals' in self._tool_inputs:
            input_string += "-L {} ".format(self._tool_inputs['TXT_intervals'][0].path)

        # set reference dbSNP db
        if 'VCF_DBSNP' in self._tool_inputs:
            input_string += "--dbsnp {} ".format(self._tool_inputs['VCF_DBSNP'][0].path)
        else:
            # set default
            self.__dbsnp_path = ToolIODb('broad_b37_dbSNP-138')
            input_string += "--dbsnp {} ".format(self.__dbsnp_path)
            logging.info("Setting dbSNP reference to default: {}".format(self.__dbsnp_path))

        if 'BAM_TUMOR' in self._tool_inputs:
            input_string += "-I:tumor {} ".format(self._tool_inputs['BAM_TUMOR'][0].path)

        if 'BAM_NORMAL' in self._tool_inputs:
            input_string += "-I:normal {} ".format(self._tool_inputs['BAM_NORMAL'][0].path)

        return input_string

    def _build_options(self, excluded_parameters=None, delimiter=' '):
        """
        Builds the options string.
        :parameter delimiter: Delimiter between option and value
        :return: Options string
        """
        if self._parameters['generate_vcf_file'].value == "False":
            excluded_parameters.append('output_vcf_file')

        options = super(Mutect1, self)._build_options(excluded_parameters=excluded_parameters, delimiter=delimiter)

        return options

    def __set_output(self):
        """
        Set the output specifications in the Camel ouptut list: 
        - plots pdf file, 
        - data csv file (optional)
        Supersedes _set_output in GATK class.
        :return: None
        """
        self._tool_outputs['TXT_CALL_STATS'] = [
            ToolIOFile(os.path.join(self._folder, self._parameters['output_callstats_file'].value))]
        if self._parameters['generate_vcf_file'].value == 'True':
            self._tool_outputs['VCF'] = [
                ToolIOFile(os.path.join(self._folder, self._parameters['output_vcf_file'].value))]

    def _check_command_output(self):
        """
        Check the result of GATK tool run
        :return: None
        """
        if len(self.stdout.split('\n')) > 1:
            if not re.match('Exit status: 0', self.stdout.split('\n')[-2].rstrip()):
                raise ToolExecutionError("Mutect1 fails to run, message: \n{}".format(self.stdout))
