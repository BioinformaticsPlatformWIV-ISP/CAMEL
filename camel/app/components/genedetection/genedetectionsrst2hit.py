import json
from typing import Optional, List

import re

from camel.app.components.genedetection.genedetectionhit import GeneDetectionHit
from camel.app.components.genedetection.genedetectionutils import GeneDetectionUtils
from camel.app.components.genedetection.mapping import Mapping
from camel.app.components.html.htmltablecell import HtmlTableCell


class GeneDetectionSRST2Hit(GeneDetectionHit):
    """
    Gene detection hit detected by SRST2.
    """

    _TABLE_COLUMNS = ['Locus', 'Length', 'Coverage', 'Mismatches', 'Uncertainty', 'Depth', 'Accession']
    _HTML_COLUMNS = _TABLE_COLUMNS

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
        self.accession = accession
        self._mismatches = mismatches if mismatches != '' else '0'
        self._uncertainty = uncertainty if uncertainty != '' else '-'
        self._depth = depth
        self._coverage = coverage
        self._length = length
        self._subject = subject
        self._extra_column_value = None
        self._extra_column_name = None

    @property
    def subject(self):
        """
        Returns the subject name.
        :return: Subject
        """
        return self._subject

    @staticmethod
    def create_from_srst2_output_line(line: str, mapping: Mapping, extra_column_value: Optional[str]):
        """
        Creates a hit object from a line in the output of SRST2.
        :param line: SRST2 output line
        :param mapping: Mapping to original sequences
        :param extra_column_value: Additional column for the hit
        :return: Hit object
        """
        parts = line.split('\t')
        full_header = mapping.get(parts[3])
        m = re.match('^(.*) ({.*})$', full_header)
        allele_full = m.group(1)
        metadata = json.loads(m.group(2))
        hit = GeneDetectionSRST2Hit(allele_full, metadata['allele'], parts[6], parts[7], float(parts[5]),
                                    float(parts[4]), int(parts[9]), metadata.get('accession', '-'))
        if extra_column_value is not None:
            name, key = GeneDetectionUtils.parse_extra_column_param(extra_column_value)
            hit.set_extra_column(name, metadata[key])
        return hit

    def to_table_row(self):
        """
        Returns the hit as a table row.
        :return: Table row
        """
        row_data = [
            self.locus,
            str(self._length),
            '{:.2f}'.format(self._coverage),
            self._mismatches,
            self._uncertainty,
            '{:.2f}'.format(self._depth),
            self.accession]
        if self._extra_column_value is not None:
            row_data.insert(-1, self._extra_column_value)
        return '\t'.join(row_data)

    def to_html_row(self, base_dir=None, sub_dir=None):
        """
        Returns the hit as a HTML row.
        :param base_dir: Base directory to store report
        :param sub_dir: Specific subdirectory of the base directory to store report files
        :return: Table row
        """
        html_data = [
            self.locus,
            self._length,
            '{:.2f}'.format(self._coverage),
            self._mismatches,
            self._uncertainty,
            '{:.2f}'.format(self._depth)]
        if self._extra_column_value is not None:
            html_data.append(self._extra_column_value)
        return [HtmlTableCell(v, self.color) for v in html_data] + [self.get_accession_cell()]

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

    @staticmethod
    def get_column_names_html(extra_column_name: Optional[str]=None) -> List[str]:
        """
        Returns the column names for the HTML output.
        :param extra_column_name: Extra column name (None if there is None)
        :return: List of column names
        """
        if extra_column_name is None:
            return GeneDetectionSRST2Hit._HTML_COLUMNS
        columns = GeneDetectionSRST2Hit._HTML_COLUMNS.copy()
        columns.insert(-1, extra_column_name)
        return columns

    @property
    def column_names_html(self):
        """
        Returns the HTML column names.
        :return: HTML column names
        """
        return GeneDetectionSRST2Hit.get_column_names_html(self._extra_column_name)

    def get_table_column_names(self):
        """
        Returns the table column names.
        :return: Table column names
        """
        if self._extra_column_name is None:
            return GeneDetectionSRST2Hit._TABLE_COLUMNS
        columns = GeneDetectionSRST2Hit._TABLE_COLUMNS.copy()
        columns.insert(-1, self._extra_column_name)
        return columns

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
