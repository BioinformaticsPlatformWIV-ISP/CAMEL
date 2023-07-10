from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.tools.blast.blast import Blast


class RpsBlast(Blast):
    """
    Protein domain search
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('rpsblast', '2.13.0', camel)

    def _check_input(self) -> None:
        """
        Checks whether the required input files are specified.
        :return: None
        """
        if 'DB_BLAST' not in self._tool_inputs:
            raise InvalidInputSpecificationError('RPS BLAST Database (key: DB_BLAST) is required!')
        super(RpsBlast, self)._check_input()
