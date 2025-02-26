from camel.app.camel import Camel
from camel.app.tools.medaka.medaka import Medaka


class MedakaVcf(Medaka):

    """
    Class for Medaka vcf function.

    Runs the medaka vcf/variant algorithm and outputs a VCF file.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes Medaka vcf.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('medaka vcf', '2.0.0', camel)

        self._required_inputs = ['HDF', 'FASTA']
        self._output_type = 'VCF'

    def _set_input(self) -> None:
        """
        Sets the input specifications and the input string.
        :return: None
        """
        super()._set_input()
        hdf_file = self._tool_inputs['HDF'][0].path
        fasta_file = self._tool_inputs['FASTA'][0].path
        self._input_string = f'{hdf_file} {fasta_file}'
