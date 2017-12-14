from app.components.genedetection.genedetectionhit import GeneDetectionHit
from app.components.html.htmltablecell import HtmlTableCell


class SRST2Hit(GeneDetectionHit):
    """
    Sequence tying hit detected by SRST2.
    """

    _TABLE_COLUMNS = ['Locus', 'Length', 'Coverage', 'Mismatches', 'Uncertainty', 'Depth', 'Accession']

    def __init__(self, subject, locus, mismatches, uncertainty, depth, coverage, length, accession):
        """
        Initializes the hit.
        :param subject: Full name of the subject sequence
        :param locus: Locus
        :param mismatches: Mismatches between reads and allele
        :param uncertainty: Uncertainty between reads and allele
        :param depth: mean read depth across allele
        :param coverage: Indicates the % of gene that was covered
        :param length: Locus length
        :param accession: Accession number
        """
        super().__init__(locus)
        self._mismatches = mismatches if mismatches != '' else '0'
        self._uncertainty = uncertainty if uncertainty != '' else '-'
        self._depth = depth
        self._coverage = coverage
        self._length = length
        self._accession = accession
        self._subject = subject

    @property
    def subject(self):
        """
        Returns the subject name.
        :return: Subject
        """
        return self._subject

    @staticmethod
    def create_from_srst2_output_line(line, mapping, metadata):
        """
        Creates a hit object from a line in the output of SRST2.
        :param line: SRST2 output line
        :param mapping: Mapping to original sequences
        :param metadata: Database metadata
        :return: Hit object
        """
        parts = line.split('\t')
        allele_full = mapping[parts[3]].split(' ')[0]
        allele_metadata = metadata[allele_full]
        return SRST2Hit(allele_full, allele_metadata['allele'], parts[6], parts[7], parts[5], float(parts[4]),
                        int(parts[9]), allele_metadata.get('accession', '-'))

    def to_table_row(self):
        """
        Returns the hit as a table row.
        :return: Table row
        """
        return "\t".join([
            self.locus,
            str(self._length),
            str(self._coverage),
            self._mismatches,
            self._uncertainty,
            str(self._depth),
            self._accession])

    def to_html_row(self, base_dir=None, sub_dir=None):
        """
        Returns the hit as a HTML row.
        :param base_dir: Base directory to store report
        :param sub_dir: Specific subdirectory of the base directory to store report files
        :return: Table row
        """
        return [HtmlTableCell(t, self.color) for t in [
            self.locus,
            self._length,
            self._coverage,
            self._mismatches,
            self._uncertainty,
            self._depth]] + [self.get_accession_cell()]

    @property
    def color(self):
        """
        Returns the hit color.
        Green: No mismatches and completely covered
        Light green: Completely covered with mismatches
        Grey: Not completely covered and mismatches
        :return: Color
        """
        if self._mismatches == '0' and self._coverage == 100.0:
            return 'green'
        elif self._coverage == 100.0:
            return 'lightgreen'
        else:
            return 'grey'

    def get_html_column_names(self):
        """
        Returns the column names for the HTML output.
        :return: Column names
        """
        return SRST2Hit._TABLE_COLUMNS

    def get_table_column_names(self):
        """
        Returns the column names for the tabular output.
        :return: Column names
        """
        return SRST2Hit._TABLE_COLUMNS
