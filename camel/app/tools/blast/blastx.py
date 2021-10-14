from camel.app.camel import Camel
from camel.app.tools.blast.blast import Blast


class Blastx(Blast):
    """
    Nucleotide - protein BLAST.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize tool.
        :return: None
        """
        super().__init__('blastx', '2.6.0', camel)
