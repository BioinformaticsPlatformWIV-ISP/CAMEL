import os

from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError


class Oncotator(Tool):
    """
    ===========
    Oncotator: anotates vcf with cancer-centered databases. 
    ===========

    Required inputs:
    ----------------
    "IN_FILE":          ToolIOFile object. {TCGAMAF,SEG_FILE,VCF,MAFLITE}. See descriptions on https://gatkforums.broadinstitute.org/gatk/discussion/4177/what-input-files-can-i-annotate-with-oncotator

    Output:
    -------
    "OUT_FILE":         ToolIOFile object. {TCGAMAF,VCF,SIMPLE_TSV,TCGAVCF,SIMPLE_BED,GENE_LIST}

    Mandatory parameters:
    ---------------------
    - db_dir        Location of databases to use for annotation
                    default: None
    - output_file   Output filename.
                    default:  oncotator.tsv

    Other parameters:
    -----------------
    - canonical_tx_file: 
                    File to use to supersede specific variants known to be clinically relevant
                    Default: None
    """

    def __init__(self, camel):
        """
        Initialize Mutect1 tool.
        :param camel: Camel instance
        :return: None
        """
        super(Oncotator, self).__init__('Oncotator', '1.9.9.0', camel)
        self._required_inputs = ['VCF']

    def _execute_tool(self):
        """
        Runs tool
        :return: None
        """
        self.__build_command()
        self._execute_command()
        self.__set_output()

    def _check_parameters(self):
        """
        Checks that parameters are valid.
        :return: None 
        """
        super(Oncotator, self)._check_parameters()

    def __build_command(self):
        """
        Build the command to run the tool.
        :return: 
        """
        input_string = " {}".format(self._tool_inputs['IN_FILE'][0].path)
        output_string = " {}".format(self._parameters['output_file_name'].value)
        options_string = " ".join(self._build_options(excluded_parameters=["output_file_name"]))

        self._command.command = "{command} {options_string} {input_string} {output_string} hg19".format(command=self._tool_command, options_string=options_string, input_string=input_string, output_string=output_string)

    def __set_output(self):
        """
        Set the output specifications in the Camel ouptut list: 
        - output file
        :return: None
        """
        self._tool_outputs['OUT_FILE'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output_file_name'].value))]

    def _check_command_output(self):
        """
        Check the result of tool run
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Oncotator failed to run with message: \n{}".format(self.stderr))
