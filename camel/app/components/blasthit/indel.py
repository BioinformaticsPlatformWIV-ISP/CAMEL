from __future__ import annotations

from camel.app.components.blasthit.blastnhit import BlastnHit


class Indel(object):

    """
    Class to represent a indel

    IMPORTANT NOTES:
     - 1-based position coordinate (first base is 1, instead of 0 (0-based))
     - the start of a indel is the position of the base proceeding it
     - deletion map to a section of REF bases between (start + 1, start + length)
     - insertion map to only 1 REF base where it starts (after this base)
     - length is the effective length of an indel = number of bases deleted or inserted

    Coordinate reference:
     - http://genomewiki.ucsc.edu/index.php/Coordinate_Transforms
    """

    def __init__(self, indel_type: str, start: int, length: int) -> None:
        """
        Initialize an Indel object
        :param indel_type: Indel type, either + or -
        :param start: Start position of the indel (1-based coordinate)
        :param length: Lenght of the indel
        :return: None
        """
        self.type = indel_type
        self.start = start
        if self.type == '+':
            self.length = length
        else:
            self.length = -length
        self.id = str(self)

    def __hash__(self) -> int:
        """
        Overwrite hash function, required for membership check, e.g., 'in'
        :return: hash value of self.id
        """
        return hash(self.id)

    def __eq__(self, other) -> bool:
        """
        Overwrite eq check (==) operation
        :param other: the other Indel to compare with
        """
        if isinstance(other, Indel):
            if self.start == other.start and self.length == other.length and self.type == other.type:
                return True

        return False

    def __str__(self) -> str:
        """
        Customized conversion to string
        :return: string presentation of Indel object
        """
        return f"({self.type}, {self.start}, {self.length})"

    def __repr__(self) -> str:
        """
        Customized output to work with 'format' function
        :return: string presentation of Indel object
        """
        # reference: http://stackoverflow.com/questions/1436703/difference-between-str-and-repr-in-python
        return str(self)

    def overlap(self, indel: Indel) -> bool:
        """
        Check if the given indel overlaps with self. Note that both indels are positioned on the query seq
        :param indel: an Indel object to compare with
        :return: boolean, True if overlap, False otherwise
        """
        if self.type == '-' and indel.type == '-':
            # both deletions
            # self:         ----
            # indel:      --    -----
            if indel.start + abs(indel.length) <= self.start or self.start + abs(self.length) <= indel.start:
                return False
            else:
                return True

        elif self.type == '-' and indel.type == '+':
            # self deletion, indel insertion
            # self:         ----
            # indel:    ++++    ++++
            if self.start >= indel.start + 1 or self.start + abs(self.length) <= indel.start:
                return False
            else:
                return True

        elif self.type == '+' and indel.type == '-':
            # self insertion, indel deletion
            # self:         +++
            # indel:     ---   ----
            if indel.start + abs(indel.length) <= self.start or self.start <= indel.start:
                return False
            else:
                return True

        elif self.type == '+' and indel.type == '+':
            # both insertions
            # self:         +++
            # indel:    ++++   ++
            if self.start == indel.start:
                return True
            else:
                return False

    def hit_overlap(self, hit: BlastnHit) -> bool:
        """
        Check if self overlaps with a blastn hit, note that both are positioned on the query seq
        :param hit: a blastn hit
        :return: boolean, True if overlap, False otherwise
        """
        # set hit as a deletion type of indel then use 'overlap' to check
        return self.overlap(Indel('-', hit.qstart - 1, hit.qend - hit.qstart + 1))

    def compatible(self, indel: Indel) -> bool:
        """
        Check if input indel is compatible with self
        :param indel: an Indel object to compare with
        :return: boolean, True if compatible, False otherwise
        """
        # compatible = exactly the same
        return self == indel
