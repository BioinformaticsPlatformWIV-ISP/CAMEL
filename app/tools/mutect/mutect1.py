import os
import re
import logging

from app.io.tooliofile import ToolIOFile
from app.io.tooliodb import ToolIODb
from app.tools.tool import Tool
from app.error.toolexecutionerror import ToolExecutionError
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError


class Mutect1(Tool):
    """
    Class for the tool Mutect (v1). Mutect performs variant calling for oncology-related NGS data. 
    V1 only calls snps, not indels.
    """

    def __init__(self, camel):
        """
        Initialize Mutect1 tool.
        :param camel: Camel instance
        :return: None
        """
        super(Mutect1, self).__init__('Mutect1', '1.1.7', camel)

        self._required_inputs = ['BAM_TUMOR']

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

    def __build_command(self):
        """
        Build the command to run the tool.
        Concatenates strings into the command string.
        :return: 
        """
        input_string = self.__create_input_string()
        options_string = ' '.join(self._build_options())
        self._command.command = ' '.join([self._tool_command, input_string, options_string])

    def __create_input_string(self):
        """
        Add the input specification in the input_string: 
        - DBSNP VCF file
        - reference fasta (default or superseded)
        - tumour BAM file
        - normal tissue BAM file (optional)
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

        # set reference dbSNP db
        if 'DBSNP_VCF' in self._tool_inputs:
            input_string += "--dbsnp {} ".format(self._tool_inputs['DBSNP_VCF'][0].path)
        else:
            # set default
            self.__dbsnp_path = ToolIODb('broad_b37_dbSNP-138')
            input_string += "--dbsnp {} ".format(self.__dbsnp_path)
            logging.info("Setting dbSNP reference to default: {}".format(self.__dbsnp_path))

        if 'TUMOR_BAM' in self._tool_inputs:
            input_string += "-I:tumor {} ".format(self._tool_inputs['TUMOR_BAM'][0].path)

        if 'NORMAL_BAM' in self._tool_inputs:
            input_string += "-I:normal {} ".format(self._tool_inputs['NORMAL_BAM'][0].path)

        return input_string

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
        if 'output_vcf_file' in self._parameters:
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
