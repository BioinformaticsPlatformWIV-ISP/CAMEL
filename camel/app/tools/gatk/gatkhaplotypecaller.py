import os

from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.gatk.gatk import GATK


class GATKHaplotypeCaller(GATK):

    """
    Class for GATK HaplotypeCaller function
    """

    def __init__(self, camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super(GATKHaplotypeCaller, self).__init__('gatk HaplotypeCaller', '3.7', camel)

        self._function_name = 'HaplotypeCaller'
        self._required_inputs = ['BAM', 'FASTA_REF']
        self._output_type = 'VCF'

    def _set_input(self):
        """
        Set the input specification
        :return: None
        """
        super(GATKHaplotypeCaller, self)._set_input()

        bam_file = self._tool_inputs['BAM'][0].path
        self._input_string += "-I {} ".format(bam_file)

    def _set_output(self):
        """
        Set the output specification
        :return: None
        """
        super(GATKHaplotypeCaller, self)._set_output()

        if 'bamOutput' in self._parameters:
            self._tool_outputs['BAM'] = [ToolIOFile(os.path.join(self._folder, self._parameters['bamOutput'].value))]
