from app.tools.blast.blast import Blast


class Blastx(Blast):
    """
    Nucleotide - protein BLAST.
    """

    def __init__(self, camel):
        """
        Initialize tool.
        :return: None
        """
        super(Blastx, self).__init__('blastx', '2.6.0', camel)
