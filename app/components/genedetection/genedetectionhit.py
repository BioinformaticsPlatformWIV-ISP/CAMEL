import abc
from abc import ABC

from app.components.html.htmltablecell import HtmlTableCell


class GeneDetectionHit(ABC):
    """
    This class represents a gene detection hit.
    """

    def __init__(self, locus):
        """
        Initializes the typing hit.
        :param locus: Locus
        """
        self._locus = locus
        self._accession = None

    @property
    def locus(self):
        """
        Returns the hit locus.
        :return: Locus
        """
        return self._locus

    @locus.setter
    def locus(self, locus):
        """
        Sets the locus.
        :param locus: Locus
        :return: None
        """
        self._locus = locus

    @property
    def accession(self):
        """
        Returns the accession number of the hit.
        :return: Accession number
        """
        return self._accession

    @accession.setter
    def accession(self, accession):
        """
        Sets the accession number of the hit.
        :param accession: Accession
        :return: None
        """
        self._accession = accession

    @abc.abstractmethod
    def to_table_row(self):
        """
        Returns the hit as a row in a table.
        :return: Table row
        """
        pass

    @abc.abstractmethod
    def to_html_row(self, report_section, sub_directory):
        """
        Returns the hit as a row in a HTML table.
        :param report_section: Section is passed to save the alignments
        :param sub_directory: Subdirectory to save the alignments
        :return: Table row
        """
        pass

    @abc.abstractmethod
    def get_table_column_names(self):
        """
        Returns the table column names.
        :return: Table column names
        """
        pass

    @abc.abstractmethod
    def get_html_column_names(self):
        """
        Returns the HTML column names.
        :return: HTML column names
        """
        pass

    @property
    @abc.abstractmethod
    def color(self):
        """
        Color for the hit.
        :return: Color
        """
        pass

    def get_accession_cell(self):
        """
        Returns the table cell for the accession.
        :return: Table cell.
        """
        if self._accession is None:
            return '-'
        else:
            link = 'https://www.ncbi.nlm.nih.gov/nuccore/{}'.format(self._accession)
            return HtmlTableCell(self._accession, self.color, link=link)
