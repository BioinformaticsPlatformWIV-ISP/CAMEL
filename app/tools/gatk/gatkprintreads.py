import os

from app.io.tooliofile import ToolIOFile
from app.tools.gatk.gatk import GATK
from app.io.tooliodb import ToolIODb
import logging


class GATKPrintReads(GATK):
    """
    Class for the GATK PrintReads tool.
    """

    def __init__(self, camel):
        """
        Initialize GATKPrintReads tool.
        :param camel: Camel instance
        :return: None
        """
        super(GATKPrintReads, self).__init__('gatk PrintReads', '3.7', camel)

        self._function_name = 'PrintReads'
        self._required_inputs = ['BAM', 'BQSR']


    def _set_input(self):
        """
        Set the input specification in the input_string
        :return: None
        """

        # set input BAM
        self._input_string += "-I {} ".format(self._tool_inputs['BAM'][0].path)

        # set input recalibration table
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
        Set the output specification in the output_string:
        BAM and BAI files
        Supersedes the _set_output fct in GATK class.
        :return: None
        """
        self.__bam_output_file = self._parameters['bam_output'].value
        self._tool_outputs['BAM'] = [ToolIOFile(os.path.join(self._folder, self.__bam_output_file))]
        self.__bai_output_file = self.__bam_output_file[:-1]+"i"
        self._tool_outputs['BAI'] = [ToolIOFile(os.path.join(self._folder, self.__bai_output_file))]
