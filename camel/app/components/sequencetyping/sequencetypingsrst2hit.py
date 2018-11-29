from typing import Optional

from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.components.sequencetyping.sequencetypinghit import SequenceTypingHit


class SequenceTypingSRST2Hit(SequenceTypingHit):
    """
    Sequence tying hit detected by SRST2.
    """

    _TABLE_COLUMNS = ['Locus', 'Allele', 'Mismatches', 'Uncertainty', 'Depth']

    def __init__(self, locus, allele_id, mismatches, uncertainty, depth):
        """
        Initializes the hit.
        :param locus: Locus
        :param allele_id: Allele id
        :param mismatches: Mismatches between reads and allele
        :param uncertainty: Uncertainty between reads and allele
        :param depth: mean read depth across allele
        """
        super().__init__(locus, allele_id)
        self._mismatches = mismatches
        self._uncertainty = uncertainty
        self._depth = depth

    @staticmethod
    def create_empty_hit(locus):
        """
        Creates an empty SRST2 hit.
        :param locus: Locus
        :return: SRST2 hit
        """
        return SequenceTypingSRST2Hit(locus, '-', '-', '-', '-')

    def to_table_row(self, separator: Optional[str]='\t'):
        """
        Returns the hit as a table row.
        :param separator: Separator for the table row
        :return: Table row
        """
        return separator.join([
            self.locus,
            self.allele_id,
            self._mismatches,
            self._uncertainty,
            '{:.2f}'.format(float(self._depth)) if self._depth != '-' else '-'
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
            self.get_allele_id_cell(),
            self._mismatches,
            self._uncertainty,
            '{:.2f}'.format(float(self._depth)) if self._depth != '-' else '-'
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

    def get_allele_id_cell(self):
        """
        Returns the cell containing the allele id.
        :return: HTML cell
        """
        if self.allele_id == '-':
            color = 'red'
        elif self.is_perfect_hit():
            color = 'green'
        elif self._uncertainty == "-":
            color = 'lightgreen'
        else:
            color = 'grey'
        return HtmlTableCell(self.allele_id, color, link=self.allele_page_url)

    def is_perfect_hit(self):
        """
        Returns true if this is a perfect hit.
        :return: True if perfect
        """
        return self._mismatches == "0" and self._uncertainty == "-"

    def __repr__(self) -> str:
        """
        Returns the internal representation.
        :return: Representation
        """
        return f"TypingSrst2Hit('{self.locus}', allele_id='{self.allele_id}', mismatch='{self._mismatches}')"
