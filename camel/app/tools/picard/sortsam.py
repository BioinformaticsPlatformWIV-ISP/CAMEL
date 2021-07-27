from camel.app.camel import Camel
from camel.app.tools.picard.picard import Picard


class SortSam(Picard):

    """
    Class for Picard SortSam function
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard SortSam', '2.23.3', camel)
        self._function_name = 'SortSam'
        self._extra_inputs = ["FASTA_REF"]

    def _set_input(self) -> None:
        """
        Set input specification
        :return: None
        """
        super(SortSam, self)._set_input()

        # optional
        if 'FASTA_REF' in self._tool_inputs:
            self._input_string += f'R={self._tool_inputs["FASTA_REF"][0].path} '
