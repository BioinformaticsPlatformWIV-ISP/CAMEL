from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.tools.picard.picard import Picard


class GatherBamFiles(Picard):

    """
    Class for Picard GatherBamFiles function.
    Concatenate efficiently BAM files that resulted from a scattered parallel analysis.
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard GatherBamFiles', '2.23.3', camel)
        self._function_name = 'GatherBamFiles'
        self._main_inputs = ['BAMs']

    def _set_input(self) -> None:
        """
        Checks length and sets input
        Overrides method in parent class.
        :return: None
        """
        input_files = [f.path for f in self._tool_inputs['BAMs']]

        if len(input_files) <= 1:
            raise InvalidInputSpecificationError("Picard GatherBamFiles: more than 1 input BAM file is expected")

        self._input_string += " I="
        self._input_string += " I=".join(input_files)