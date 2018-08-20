import os

from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.components.sequencetyping.sequencetypinghit import SequenceTypingHit


class SequenceTypingBlastHit(SequenceTypingHit):
    """
    Sequence tying hit detected by blast.
    """

    _TABLE_COLUMNS = ['Locus', 'Allele', '% Identity', 'HSP/Locus length', 'Type']
    _HTML_COLUMNS = _TABLE_COLUMNS + ['Alignment']

    def __init__(self, locus, allele_id, type_, subject, pident, slen, sseq, qseqid, qstart, qend):
        """
        Initializes the hit.
        :param locus: Locus
        :param allele_id: Allele id
        :param type_: Locus type ('DNA', 'peptide')
        :param subject: Subject
        :param pident: Percent identity
        :param slen: Subject length
        :param sseq: Aligned sequence of the subject
        :param qseqid: Query sequence id
        :param qstart: Start of the alignment in the query
        :param qend: End of the alignment in the query
        """
        super().__init__(locus, allele_id)
        self._subject = subject
        self._pident = float(pident) if pident != '-' else None
        self._type = type_
        self._slen = slen
        self._sseq = sseq
        self._qsedid = qseqid
        self._qstart = qstart
        self._qend = qend
        self._alignment_path = None

    @staticmethod
    def create_from_dict(input_dict, type_):
        """
        Creates a hit object from a dictionary containing the blast output.
        Allele id is set to None, it extracted afterwards only for the best hits.
        :param input_dict: Input dictionary
        :param type_: Locus type
        :return: Hit object
        """
        try:
            return SequenceTypingBlastHit(None, None, type_, input_dict['sseqid'], float(input_dict['pident']),
                                          input_dict['slen'], input_dict['sseq'], input_dict['qseqid'],
                                          input_dict['qstart'], input_dict['qend'])
        except KeyError as err:
            raise ValueError("Cannot create hit from dictionary {} missing - {!r}".format(err, input_dict))

    def to_table_row(self):
        """
        Converts the hit into a table row.
        :return: Table row
        """
        return '\t'.join([
            self.locus,
            self.allele_id,
            '{:.2f}'.format(self.percent_identity) if self.percent_identity is not None else '-',
            self.length_statistic,
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
            '{:.2f}'.format(self.percent_identity) if self.percent_identity is not None else '-',
            self.length_statistic,
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
    def generate_empty_hit(locus, type_):
        """
        Returns an empty hit.
        :param locus: Locus
        :param type_: Locus type
        :return: None
        """
        return SequenceTypingBlastHit(locus, '-', type_, '-', '-', '-', '-', '-', '-', '-')

    @staticmethod
    def generate_multi_hit(locus, type_):
        """
        Returns a multi hit.
        :param locus: Locus
        :param type_: Locus type
        :return: None
        """
        return SequenceTypingBlastHit(locus, '?', type_, '-', '-', '-', '-', '-', '-', '-')

    @property
    def subject(self):
        """
        Returns the subject (locus + allele id).
        :return: Subject
        """
        return self._subject

    @property
    def query(self):
        """
        Returns the query.
        :return: Query
        """
        return self._qsedid

    @property
    def subject_length(self):
        """
        Returns the subject length.
        :return: Subject length
        """
        return self._slen

    @property
    def alignment_length(self):
        """
        Returns the alignment length.
        :return: Alignment length
        """
        return len(self._sseq)

    @property
    def percent_identity(self):
        """
        Returns the percent identity.
        :return: Percent identity
        """
        return self._pident

    @property
    def subject_coverage(self):
        """
        Returns the fraction of the subject that is covered by the alignment.
        :return: % subject covered
        """
        return 100.0 * float(self.alignment_length) / self._slen

    @property
    def length_statistic(self):
        """
        Returns the subject coverage in the format: {bases_covered}/{subject_length}.
        :return: Length statistic
        """
        if self._slen == '-':
            return '-'
        return '{}/{}'.format(self.alignment_length, self._slen)

    @property
    def gaps(self):
        """
        Returns the number of gaps in the alignment.
        :return: Number of gaps
        """
        return self._sseq.count('-')

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

    def is_perfect_hit(self):
        """
        Returns true if the hit is perfect (100% identity over complete length)
        :return: True if perfect
        """
        return (self.percent_identity == 100.0) and (self._slen == self.alignment_length)

    @property
    def color(self):
        """
        Returns the color for this hit.
        Green: Perfect hit
        Light green: Full length hit with one or more mismatches
        Grey: Non-full length hit
        Red: No-hit
        :return: Color
        """
        if self.is_perfect_hit():
            return 'green'
        elif self._slen == self.alignment_length:
            return 'lightgreen'
        elif self.percent_identity != '-':
            return 'grey'
        elif self.allele_id == '?':
            return 'yellow'
        else:
            return 'red'
