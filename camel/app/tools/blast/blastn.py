from camel.app.camel import Camel
from camel.app.tools.blast.blast import Blast


class Blastn(Blast):
    """
    Nucleotide - nucleotide BLAST.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize tool.
        :return: None
        """
        super().__init__('blastn', '2.14.0', camel)
