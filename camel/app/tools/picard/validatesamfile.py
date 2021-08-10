from camel.app.camel import Camel
from camel.app.tools.picard.picard import Picard


class ValidateSamFile(Picard):
    """
    Class for Picard ValidateSamFile function
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard ValidateSamFile', '2.23.3', camel)
        self._output_type = 'TXT_report'
