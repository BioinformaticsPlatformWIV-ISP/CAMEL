class SnpPosition(object):
    """
    Class that represent a position in the genome where a SNP occurred.
    """

    def __init__(self, contig, position, reference_base):
        """
        Initializes a SNP position.
        :param contig: Contig
        :param position: Position in the contig
        """
        self.contig = contig
        self.position = position
        self.reference_base = reference_base

    def __eq__(self, other):
        """
        Overrides equality.
        :param other: Other position
        :return: True if the positions are the same
        """
        if (self.position == other.position) and (self.contig == other.contig):
            return True
        return False

    def __hash__(self):
        """
        Hash function, this ensures that SNP's at the same location get the same key.
        :return: Hash
        """
        return hash((self.contig, self.position))

    def __lt__(self, other: 'SnpPosition') -> bool:
        """
        Compares two SNP positions.
        :param other: Other SNP position
        :return: True if the other is larger
        """
        if self.contig == other.contig:
            return self.position < other.position
        else:
            return self.contig > other.contig
