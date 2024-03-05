from camel.app.camel import Camel
from camel.app.tools.blast.blast import Blast


class Blastp(Blast):
    """
    BlastP compares a protein query to a protein database.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize tool.
        :return: None
        """
        super().__init__('blastp', '2.14.0', camel)
