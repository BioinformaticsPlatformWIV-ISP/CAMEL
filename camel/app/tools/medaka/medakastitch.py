from camel.app.camel import Camel
from camel.app.tools.medaka.medaka import Medaka


class MedakaStitch(Medaka):

    """
    Class for Medaka stitch function.

    Medaka stitch reads the output of Medaka consensus and outputs a consensus fasta sequence.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes Medaka stitch.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('medaka stitch', '1.7.3', camel)

        self._required_inputs = ['HDF', 'FASTA']
        self._output_type = 'FASTA'

    def _set_input(self) -> None:
        """
        Sets the input specification and the input string.
        :return: None
        """
        super()._set_input()
        hdf_file = self._tool_inputs['HDF'][0].path
        fasta_file = self._tool_inputs['FASTA'][0].path
        self._input_string = f'{hdf_file} {fasta_file}'
