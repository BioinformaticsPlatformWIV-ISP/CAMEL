import os

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.picard.picard import Picard


class MergeVCFs(Picard):
    """
    Class for picard MergeVCFs function
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard MergeVCFs', '2.23.3', camel)

        self._function_name = 'MergeVCFs'
        self._main_inputs = ['VCF']

    def _set_input(self) -> None:
        """
        Set the input specification. This method handles on or more VCF files
        :return: None
        """
        self._input_string += "".join(["I=", " I=".join([vcf.path for vcf in self._tool_inputs["VCF"]])])

    def _set_output(self) -> None:
        """
        Set the output specification, this default function handles one VCF file as output
        :return: None
        """
        self._tool_outputs['VCF'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output'].value))]
