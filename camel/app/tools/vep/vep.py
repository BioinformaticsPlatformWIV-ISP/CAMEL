import os
import re

from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError


class Vep(Tool):
    """
    ===========
    
    ===========
    
    Required inputs:
    ----------------
    "VCF":              ToolIOFile object. VCF or TXT (whitespace seperated) file with variants.

    Optional input:
    ---------------
    "BAM_NORMAL":       ToolIOFile object. BAM file with normal data for tumor-normal matching.
    "FASTA_REF":        ToolIOFile object. FASTA file containing the reference genome.
    "VCF_DBSNP":        ToolIOFile object. DbSNP reference vcf file location.
    "TXT_intervals":    ToolIOFile object. Intervals list to restrict search by GATK. Accelerates analysis. Bed or GATK intervals list 

    Output:
    -------
    "VCF":              ToolIOFile object. VCF file with annotations. Can also be txt file.
    "HTML":             ToolIOFile object. HTML file with summary of variants.

    Mandatory parameters:
    ---------------------
    - output_file
                    default value:  vep_annotated.vcf
    
    Other parameters:
    -----------------
    - output_in_vcf    output in vcf format instead of default text format
    """

    def __init__(self, camel):
        """
        Initialize Mutect1 tool.
        :param camel: Camel instance
        :return: None
        """
        super(Vep, self).__init__('Vep', '93', camel)
        self._required_inputs = ['VCF', 'DB_PATH']

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
        super(Vep, self)._check_input()

        for input_key in self._required_inputs:
            if input_key not in self._tool_inputs:
                raise InvalidInputSpecificationError(
                    'Vep required {} input is missing in tool inputs!'.format(input_key))

    def _check_parameters(self):
        """
        Checks that parameters are valid.
        :return: None 
        """
        super(Vep, self)._check_parameters()

    def __build_command(self):
        """
        Build the command to run the tool.
        By default, use cache provided in DB_PATH and run offline.
        Concatenates strings into the command string.
        :return: 
        """
        input_string = " -i {}".format(self._tool_inputs['VCF'][0].path)

        options_string = "--cache --offline "
        options_string += " ".join(self._build_options())
        options_string += " --dir_cache {} ".format(self._tool_inputs['DB_PATH'][0].path)

        self._command.command = ' '.join([self._tool_command, input_string, options_string])

    def __set_output(self):
        """
        Set the output specifications in the Camel ouptut list: 
        - html file
        - vcf file
        :return: None
        """
        self._tool_outputs['VCF'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output_file'].value))]

        html_file_name = "{}_summary.html".format(os.path.join(self._folder, self._parameters['output_file'].value))
        self._tool_outputs['HTML'] = [ToolIOFile(html_file_name)]

    def _check_command_output(self):
        """
        Check the result of Vep tool run
        :return: None
        """
        if not self.stdout == "" and "ERROR" in self.stdout:
            raise ToolExecutionError("Vep fails to run, message: \n{}".format(self.stdout))
