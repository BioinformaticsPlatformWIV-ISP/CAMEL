import os

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.gatk4.gatk4 import GATK4


class GATK4HaplotypeCaller(GATK4):

    """
    Class for GATK4 HaplotypeCaller function
    """

    def __init__(self, camel: Camel):
        """
        Initialize the GATK4 HaplotypeCaller
        :param camel: Camel instance
        :return: None
        """
        super().__init__('gatk4 HaplotypeCaller', '4.1.9.0', camel)

        self._required_inputs = ['BAM', 'FASTA_REF']
        self._output_type = 'VCF'

    def _set_input(self) -> None:
        """
        Set the input specification
        :return: None
        """
        super(GATK4HaplotypeCaller, self)._set_input()

        bam_file = self._tool_inputs['BAM'][0].path
        self._input_string += f"--input {bam_file}"

    def _set_output(self) -> None:
        """
        Set the output specification
        :return: None
        """
        super(GATK4HaplotypeCaller, self)._set_output()

        if 'bam-output' in self._parameters:
            self._tool_outputs['BAM'] = [ToolIOFile(os.path.join(self._folder, self._parameters['bam-output'].value))]
