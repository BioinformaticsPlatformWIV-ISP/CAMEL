from typing import Optional

from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.components.sequencetyping.sequencetypinghit import SequenceTypingHit


class SequenceTypingKMAHit(SequenceTypingHit):
    """
    Sequence tying hit detected by SRST2.
    """

    _TABLE_COLUMNS = ['Locus', 'Allele', 'Length', '% Identity', '% Coverage', 'Depth']

    def __init__(self, locus: str, allele_id: str, length: int, p_ident: float, p_cov: float, depth: float,
                 score: int) -> None:
        """
        Initializes the hit.
        :param locus: Locus
        :param allele_id: Allele id
        :param length: (Subject) length
        :param p_ident: (Subject) percent identity
        :param p_cov: (Subject) percent coverage
        :param depth: mean read depth across allele
        :param score: k-mer score
        """
        super().__init__(locus, allele_id)
        self._length = length
        self._p_ident = p_ident
        self._p_cov = p_cov
        self._depth = depth
        self._score = score

    @staticmethod
    def create_empty_hit(locus):
        """
        Creates an empty KMA hit.
        :param locus: Locus
        :return: KMA hit
        """
        return SequenceTypingKMAHit(locus, '-', 0, 0, 0, 0, 0)

    def to_table_row(self, separator: Optional[str] = '\t'):
        """
        Returns the hit as a table row.
        :param separator: Separator for the table row
        :return: Table row
        """
        return separator.join([
            self.locus,
            self.allele_id,
            str(self._length),
            '{:.2f}'.format(float(self._p_ident)),
            '{:.2f}'.format(float(self._p_cov)),
            '{:.2f}'.format(float(self._depth))
        ])

    def to_html_row(self, base_dir=None, sub_dir=None):
        """
        Returns the hit as a HTML row.
        :param base_dir: Base directory to store report
        :param sub_dir: Specific subdirectory of the base directory to store report files
        :return: Table row
        """
        return [
            self.locus,
            HtmlTableCell(self.allele_id, self.color, link=self.allele_page_url),
            str(self._length),
            '{:.2f}'.format(float(self._p_ident)),
            '{:.2f}'.format(float(self._p_cov)),
            '{:.2f}'.format(float(self._depth)),
        ]

    def get_table_column_names(self):
        """
        Returns the table column names.
        :return: Table column names
        """
        return self._TABLE_COLUMNS

    def get_html_column_names(self):
        """
        Returns the HTML column names.
        :return: HTML column names
        """
        return self._TABLE_COLUMNS

    @property
    def color(self) -> str:
        """
        Returns the color for this hit.
        :return: Color
        """
        if self.allele_id == '-':
            return 'red'
        elif self.is_perfect_hit():
            return 'green'
        elif self.is_full_length():
            return 'lightgreen'
        else:
            return 'grey'

    def is_perfect_hit(self) -> bool:
        """
        Returns true if this is a perfect hit.
        :return: True if perfect
        """
        return self._p_ident == 100.0 and self.is_full_length()

    def is_full_length(self) -> bool:
        """
        Returns true if this is a full length hit.
        :return: True if full length
        """
        return self._p_cov == 100.0

    def __repr__(self) -> str:
        """
        Returns the internal representation.
        :return: Representation
        """
        return f"TypingKMAHit('{self.locus}', allele_id='{self.allele_id}')"

    @property
    def score(self) -> int:
        """
        Returns the k-mer score.
        :return: Score
        """
        return self._score

    @property
    def percent_identity(self) -> float:
        """
        Returns the percent identity of the hit.
        :return: Percent identity
        """
        return self._p_ident

    @property
    def percent_covered(self) -> float:
        """
        Returns the percentage of the subject sequence that is covered.
        :return:
        """
        return self._p_cov
