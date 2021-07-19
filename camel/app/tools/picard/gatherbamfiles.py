from camel.app.camel import Camel
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

        self._required_inputs = ['BAMs']
        self._supported_inputs = ['BAMs']

    def _check_input(self) -> None:
        """
        Check and set the input.
        Overrides method in parent class.
        :return: None
        """
        super(Picard, self)._check_input()

        self._set_input()

    def _set_input(self) -> None:
        """
        Set the input specification in the input_string
        Overrides method in parent class.
        :return: None
        """

        # set input reports
        self._input_string = ' I='
        self._input_string += ' I='.join(f.path for f in self._tool_inputs['BAMs'])