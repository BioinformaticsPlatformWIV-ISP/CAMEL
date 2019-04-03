from camel.app.tools.picard.picard import Picard


class SortSam(Picard):

    """
    Class for Picard SortSam function
    """

    def __init__(self, camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard SortSam', '2.8.3', camel)
        self._function_name = 'SortSam'
