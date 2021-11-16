from camel.app.camel import Camel
from camel.app.tools.picard.picard import Picard


class CollectQualityYieldMetrics(Picard):
    """
    Class for picard CollectQualityYieldMetrics function to calculate a set of metrics used to describe the
    general quality of a BAM file
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard CollectQualityYieldMetrics', '2.23.3', camel)
        self._output_type = 'TXT'
