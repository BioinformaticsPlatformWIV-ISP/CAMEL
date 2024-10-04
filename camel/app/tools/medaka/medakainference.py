from camel.app.camel import Camel
from camel.app.tools.medaka.medaka import Medaka


class MedakaInference(Medaka):

    """
    Class for Medaka inference function.

    Runs the medaka inference algorithm and outputs a HDF file in preparation for building a consensus sequence.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes Medaka inference.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('medaka inference', '2.0.0', camel)

        self._required_inputs = ['BAM']
        self._output_type = 'HDF'

    def _set_input(self) -> None:
        """
        Sets the input specifications and the input string.
        :return: None
        """
        super()._set_input()
        self._input_string = str(self._tool_inputs['BAM'][0].path)
