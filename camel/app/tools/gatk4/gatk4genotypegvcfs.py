from camel.app.camel import Camel
from camel.app.tools.gatk4.gatk4 import GATK4


class GATK4GenotypeGVCFs(GATK4):

    """
    Class for GATK GenotypeGVCFs function
    """

    def __init__(self, camel: Camel):
        """
        Initialize the GATK4GenotypeGVCFs tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('gatk4 GenotypeGVCFs', '4.1.9.0', camel)

        self._required_inputs = ['gVCF', 'FASTA_REF']
        self._output_type = 'VCF_MultipleSample'

    def _set_input(self) -> None:
        """
        Set the input specification
        :return: None
        """
        super(GATK4GenotypeGVCFs, self)._set_input()

        for f in self._tool_inputs['gVCF']:
            self._input_string += f"--variant {f.path} "
