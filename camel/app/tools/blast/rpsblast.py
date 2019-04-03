from camel.app.tools.blast.blast import Blast
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError


class RpsBlast(Blast):
    """
    Protein domain search
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('rpsblast', '2.6.0', camel)

    def _check_input(self):
        """
        Checks whether the required input files are specified.
        :return: None
        """
        if 'DB_BLAST' not in self._tool_inputs:
            raise InvalidInputSpecificationError('RPS BLAST Database (key: DB_BLAST) is required!')
        super(RpsBlast, self)._check_input()
