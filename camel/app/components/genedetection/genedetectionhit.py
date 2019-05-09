from typing import Any, List, Union

import abc
from abc import ABC

from camel.app.components.genedetection.genedetectionutils import GeneDetectionUtils
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell


class GeneDetectionHit(ABC):
    """
    This class represents a gene detection hit.
    """

    def __init__(self, locus: str) -> None:
        """
        Initializes the typing hit.
        :param locus: Locus
        """
        self._locus = locus
        self._accession = None

    @property
    def locus(self) -> str:
        """
        Returns the hit locus.
        :return: Locus
        """
        return self._locus

    @locus.setter
    def locus(self, locus: str) -> None:
        """
        Sets the locus.
        :param locus: Locus
        :return: None
        """
        self._locus = locus

    @property
    def accession(self) -> str:
        """
        Returns the accession number of the hit.
        :return: Accession number
        """
        return self._accession

    @accession.setter
    def accession(self, accession: str) -> None:
        """
        Sets the accession number of the hit.
        :param accession: Accession
        :return: None
        """
        self._accession = str(accession)

    @abc.abstractmethod
    def to_table_row(self) -> str:
        """
        Returns the hit as a row in a table.
        :return: Table row
        """
        pass

    @abc.abstractmethod
    def to_html_row(self, report_section: HtmlReportSection, sub_directory: str, colored: bool = True) -> List[Any]:
        """
        Returns the hit as a row in a HTML table.
        :param report_section: Section is passed to save the alignments
        :param sub_directory: Subdirectory to save the alignments
        :param colored: If True, the row is colored
        :return: Table row
        """
        pass

    @property
    @abc.abstractmethod
    def color(self) -> str:
        """
        Color for the hit.
        :return: Color
        """
        pass

    def get_accession_cell(self) -> Union[HtmlTableCell, str]:
        """
        Returns the table cell for the accession.
        :return: Table cell.
        """
        if self._accession is None:
            return '-'
        elif GeneDetectionUtils.is_ncbi_accession(self.accession):
            link = 'https://www.ncbi.nlm.nih.gov/nuccore/{}'.format(self._accession)
            return HtmlTableCell(self._accession, self.color, link=link)
        else:
            return HtmlTableCell(self._accession, self.color)
