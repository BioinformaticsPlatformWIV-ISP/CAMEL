from typing import Optional

import os

from camel.app.components.blast.blasthitstatistics import BlastHitStatistics
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.components.sequencetyping.sequencetypinghit import SequenceTypingHit


class SequenceTypingBlastHit(SequenceTypingHit):
    """
    Sequence tying hit detected by blast.
    """

    _TABLE_COLUMNS = ['Locus', 'Allele', '% Identity', 'HSP/Locus length', 'Type']
    _HTML_COLUMNS = _TABLE_COLUMNS + ['Alignment']
    SYMBOL_MULTI_HIT = '?'
    SYMBOL_NO_HIT = '-'

    def __init__(self, locus: Optional[str], allele_id: Optional[str], type_,
                 blast_stats: Optional[BlastHitStatistics]) -> None:
        """
        Initializes the hit.
        :param locus: Locus
        :param allele_id: Allele id
        :param type_: Locus type ('DNA', 'peptide')
        :param blast_stats: Blast hit statistics
        """
        super().__init__(locus, allele_id)
        self._type = type_
        self._blast_stats = blast_stats
        self._alignment_path = None

    def to_table_row(self, separator: str = '\t'):
        """
        Converts the hit into a table row.
        :param separator: Separator
        :return: Table row
        """
        return separator.join([
            self.locus,
            self.allele_id,
            '{:.2f}'.format(self._blast_stats.percent_identity) if self.blast_stats else '-',
            self.blast_stats.length_statistic if self.blast_stats else '-',
            self._type])

    def to_html_row(self, report_section, sub_dir=None):
        """
        Converts the hit into a HTML table row
        :param report_section: Section is passed to save the alignments
        :param sub_dir: Specific subdirectory of the base directory to store report files
        :return: HTML row elements
        """
        if self._alignment_path is None:
            alignment_cell = '-'
        else:
            relative_path = os.path.join(sub_dir, 'alignments', os.path.basename(self._alignment_path))
            report_section.add_file(self._alignment_path, relative_path)
            alignment_cell = HtmlTableCell('view', link=relative_path)
        return [
            self.locus,
            HtmlTableCell(self.allele_id, self.color, link=self.allele_page_url),
            '{:.2f}'.format(self.blast_stats.percent_identity) if self.blast_stats else '-',
            self.blast_stats.length_statistic if self.blast_stats else '-',
            self._type,
            alignment_cell]

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
        return self._HTML_COLUMNS

    @staticmethod
    def generate_empty_hit(locus: str, type_: str) -> 'SequenceTypingBlastHit':
        """
        Returns an empty hit.
        :param locus: Locus
        :param type_: Locus type
        :return: None
        """
        return SequenceTypingBlastHit(locus, SequenceTypingBlastHit.SYMBOL_NO_HIT, type_, None)

    @staticmethod
    def generate_multi_hit(locus: str, type_: str) -> 'SequenceTypingBlastHit':
        """
        Returns a multi hit.
        :param locus: Locus
        :param type_: Locus type
        :return: None
        """
        return SequenceTypingBlastHit(locus, SequenceTypingBlastHit.SYMBOL_MULTI_HIT, type_, None)

    @property
    def blast_stats(self) -> BlastHitStatistics:
        """
        Returns the BLAST statistics object.
        :return: BLAST hit statistics
        """
        return self._blast_stats

    @property
    def alignment_path(self):
        """
        Returns the path to the alignment file.
        :return: Alignment file
        """
        return self._alignment_path

    @alignment_path.setter
    def alignment_path(self, alignment_path):
        """
        Sets the alignment path.
        :param alignment_path: Alignment path
        :return: None
        """
        self._alignment_path = alignment_path

    @property
    def color(self) -> str:
        """
        Returns the color for this hit.
        Green: Perfect hit
        Light green: Full length hit with one or more mismatches
        Grey: Non-full length hit
        Red: No-hit
        :return: Color
        """
        if self.allele_id == SequenceTypingBlastHit.SYMBOL_MULTI_HIT:
            return 'yellow'
        elif self.allele_id == SequenceTypingBlastHit.SYMBOL_NO_HIT:
            return 'red'
        elif self.is_perfect_hit():
            return 'green'
        elif self.is_full_length():
            return 'lightgreen'
        return 'grey'

    def __repr__(self) -> str:
        """
        Returns the internal representation.
        :return: Representation
        """
        return f"TypingBlastHit('{self.locus}', allele_id='{self.allele_id}')"

    def is_perfect_hit(self) -> bool:
        """
        Function to check if this is a perfect hit.
        :return: True if perfect hit, False otherwise
        """
        if self._blast_stats is None:
            return False
        return self._blast_stats.is_perfect_hit()

    def is_full_length(self) -> bool:
        """
        Function to check if this is a full length hit.
        :return: True if full length, False otherwise
        """
        if self._blast_stats is None:
            return False
        return self._blast_stats.is_full_length()
