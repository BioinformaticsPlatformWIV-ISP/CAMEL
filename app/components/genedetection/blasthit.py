import os

from app.components.genedetection.genedetectionhit import GeneDetectionHit
from app.components.html.htmltablecell import HtmlTableCell


class BlastHit(GeneDetectionHit):
    """
    Gene detection hit detected by blast.
    """

    _TABLE_COLUMNS = ['Locus', '% Identity', 'HSP/Locus length', 'Contig', 'Position in contig', 'Accession']
    _HTML_COLUMNS = _TABLE_COLUMNS + ['Alignment']

    def __init__(self, locus, subject, pident, slen, sseq, qseqid, qstart, qend, accession, alignment_path):
        """
        Initializes the hit.
        :param locus: Locus
        :param subject: Subject
        :param pident: Percent identity
        :param slen: Subject length
        :param sseq: Aligned sequence of the subject
        :param qseqid: Query sequence id
        :param qstart: Start of the alignment in the query
        :param qend: End of the alignment in the query
        :param accession: NCBI accession number
        :param alignment_path: Path to the alignment visualization file
        """
        super().__init__(locus)
        self._subject = subject
        self._pident = pident
        self._slen = slen
        self._sseq = sseq
        self._qseqid = qseqid
        self._qstart = qstart
        self._qend = qend
        self._accession = accession
        self._alignment_path = alignment_path
        self._extra_column_value = None
        self._extra_column_name = None

    @staticmethod
    def create_from_dict(input_dict):
        """
        Creates a hit object from a dictionary containing the blast output.
        :param input_dict: Input dictionary
        :return: Hit object
        """
        try:
            return BlastHit(None, input_dict['sseqid'], float(input_dict['pident']), input_dict['slen'],
                            input_dict['sseq'], input_dict['qseqid'], input_dict['qstart'], input_dict['qend'], None,
                            None)
        except KeyError as err:
            raise ValueError("Cannot create hit from dictionary {} missing - {!r}".format(err, input_dict))

    def get_table_column_names(self):
        """
        Returns the table column names.
        :return: Table column names
        """
        if self._extra_column_name is None:
            return BlastHit._TABLE_COLUMNS
        columns = BlastHit._TABLE_COLUMNS.copy()
        columns.insert(-2, self._extra_column_name)
        return columns

    def get_html_column_names(self):
        """
        Returns the HTML column names.
        :return: HTML column names
        """
        if self._extra_column_name is None:
            return BlastHit._HTML_COLUMNS
        columns = BlastHit._HTML_COLUMNS.copy()
        columns.insert(-3, self._extra_column_name)
        return columns

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
        return self._qseqid

    @property
    def query_start(self):
        """
        Returns the start position of the query in the alignment.
        :return: Query start
        """
        return self._qstart

    @property
    def query_end(self):
        """
        Returns the end position of the query in the alignment.
        :return: Query end
        """
        return self._qend

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
        return 100.0 * float(self.alignment_length) / self.subject_length

    @property
    def length_statistic(self):
        """
        Returns the subject coverage in the format: {bases_covered}/{subject_length}.
        :return: Length statistic
        """
        return '{}/{}'.format(self.alignment_length, self.subject_length) if self.subject_length != '-' else '-'

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

    @property
    def extra_column_value(self):
        """
        Returns the value of the extra column.
        :return: Extra column
        """
        return self._extra_column_value

    def set_extra_column(self, name, value):
        """
        Sets the extra column information.
        This extra column is used to contains some additional metadata associated with this hit. It is included in
        the tabular output and the HTML output. It consists of a column name and a value.
        E.g.: name - 'Protein function', value - 'Heat shock protein'
        :param name: Name of the extra column
        :param value: Value of the extra column
        :return: None
        """
        self._extra_column_value = value
        self._extra_column_name = name

    def is_perfect_hit(self):
        """
        Returns true if the hit is perfect (100% identity over complete length)
        :return: True if perfect
        """
        return (self.percent_identity == 100.0) and (self.subject_length == self.alignment_length)

    def to_table_row(self):
        """
        Converts the hit into a table row.
        :return: Table row
        """
        row_data = [
            self.locus,
            str(self.percent_identity),
            self.length_statistic,
            self.query,
            '{}..{}'.format(self.query_start, self.query_end),
            self.accession if self.accession is not None else '-']
        if self._extra_column_value is not None:
            row_data.insert(-2, self._extra_column_value)
        return '\t'.join(row_data)

    def to_html_row(self, report_section, sub_directory):
        """
        Converts the hit into a HTML table row
        :param report_section: HTML Section that will contain the hit table
        :param sub_directory: Subdirectory to save the alignments
        :return: HTML row elements
        """
        if self.alignment_path is None:
            alignment_cell = '-'
        else:
            relative_path = os.path.join(sub_directory, 'alignments', os.path.basename(self.alignment_path))
            report_section.add_file(self.alignment_path, relative_path)
            alignment_cell = HtmlTableCell('view', self.color, link=relative_path)
        html_data = [
            self.locus,
            str(self.percent_identity),
            self.length_statistic,
            self.query,
            '{}..{}'.format(self.query_start, self.query_end)]
        if self._extra_column_value is not None:
            html_data.insert(-1, self._extra_column_value)
        return [HtmlTableCell(v, self.color) for v in html_data] + [self.get_accession_cell()] + [alignment_cell]

    @property
    def color(self):
        """
        Returns the color for this hit.
        Green: Perfect hit
        Light green: Full length hit with one or more mismatches
        Grey: Non-full length hit
        :return: Color
        """
        if self.is_perfect_hit():
            return 'green'
        elif self.subject_length == self.alignment_length:
            return 'lightgreen'
        elif self.percent_identity != '-':
            return 'grey'
