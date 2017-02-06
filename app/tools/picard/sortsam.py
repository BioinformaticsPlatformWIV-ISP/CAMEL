from app.tools.picard.picard import Picard


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
        super(SortSam, self).__init__('Picard SortSam', '2.6.0', camel)
        self.function_name = 'SortSam'
