from camel.app.camel import Camel
from camel.app.tools.gatk4.gatk4 import GATK4
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError


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
        self._specific_parameters = ['gendb']

    def _check_input(self) -> None:
        """
        Check input for a tool and prepare command line parameters for input
        :return: None
        """
        if self._parameters['gendb'] is not None:
            self._required_inputs.remove('gVCF')

        super(GATK4GenotypeGVCFs, self)._check_input()

    def _set_input(self) -> None:
        """
        Set the input specification
        :return: None
        """
        super(GATK4GenotypeGVCFs, self)._set_input()

        if 'gVCF' in self._tool_inputs:
            for f in self._tool_inputs['gVCF']:
                self._input_string += f"--variant {f.path} "

        elif 'gendb' in self._parameters:
            self._input_string += f"--variant gendb:/{self._parameters['gendb'].value}"

    def _build_command(self):
        """
        Build the command to run tool.
        :return: None
        """
        if 'annotation_group' in self._parameters:
            self._option_string += ' --annotation-group '
            self._option_string += ' --annotation-group '.join(self._parameters['annotation_group'].value.split(","))

        self._option_string += ' '

        super(GATK4HaplotypeCaller, self)._build_command()