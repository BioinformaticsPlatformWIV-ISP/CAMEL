from camel.app.camel import Camel
from camel.app.tools.medaka.medaka import Medaka


class MedakaConsensus(Medaka):

    """
    Class for Medaka consensus function.

    Runs the medaka consensus algorithm and outputs a HDF file in preparation for building a consensus sequence.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes Medaka consensus.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('medaka consensus', '1.7.3', camel)

        self._required_inputs = ['BAM']
        self._output_type = 'HDF'

    def _set_input(self) -> None:
        """
        Sets the input specifications and the input string.
        :return: None
        """
        super()._set_input()
        self._input_string = str(self._tool_inputs['BAM'][0].path)
