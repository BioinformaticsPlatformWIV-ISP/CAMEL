from camel.app.camel import Camel
from camel.app.tools.picard.picard import Picard


class SamToFastq(Picard):

    """
    Class for Picard FastqToSam function
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard SamToFastq', '2.23.3', camel)
        self._function_name = 'SamToFastq'
        self._output_type = 'FASTQ'
