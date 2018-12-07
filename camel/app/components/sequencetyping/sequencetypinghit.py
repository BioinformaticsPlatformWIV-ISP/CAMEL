import abc
from typing import Optional, List, Any

from camel.app.components.html.htmlreportsection import HtmlReportSection


class SequenceTypingHit(metaclass=abc.ABCMeta):
    """
    This class represents is the base class for sequence typing hits. All hits should define the locus to which they
    belong and the allele that was detected. Abstract methods should be implemented by sub-classes so other classes
    can rely on these methods to provide the required functionality regardless of detection method.
    """

    def __init__(self, locus: str, allele_id: str):
        """
        Initializes the typing hit.
        :param locus: Locus
        :param allele_id: Allele id of the hit
        """
        self._locus = locus
        self._allele_id = allele_id
        self._allele_page_url_template = None

    @property
    def locus(self) -> str:
        """
        Returns the hit locus.
        :return: Locus
        """
        return self._locus

    @locus.setter
    def locus(self, locus: str):
        """
        Sets the locus.
        :param locus: Locus
        :return: None
        """
        self._locus = locus

    @property
    def allele_id(self) -> str:
        """
        Returns the allele id.
        :return: Allele id
        """
        return self._allele_id

    @allele_id.setter
    def allele_id(self, allele_id: str):
        """
        Sets the allele id.
        :param allele_id: Allele id
        :return: None
        """
        self._allele_id = allele_id

    def set_allele_page_url_template(self, url_template: str) -> None:
        """
        Sets the allele page URL template.
        :param url_template: : URL template
        :return: None
        """
        self._allele_page_url_template = url_template

    @property
    def allele_page_url(self) -> Optional[str]:
        """
        Returns the allele info page URL if there is one.
        :return: Allele page URL
        """
        if self._allele_page_url_template is None:
            return None
        elif self._allele_id in ('?', '-'):
            return None
        return self._allele_page_url_template.format(allele_id=self.allele_id)

    @abc.abstractmethod
    def to_table_row(self, separator: Optional[str]='\t') -> List[Any]:
        """
        Returns the hit as a row in a table.
        :param separator: Separator
        :return: Table row
        """
        pass

    @abc.abstractmethod
    def to_html_row(self, report_section: HtmlReportSection, sub_dir: str=None) -> List[Any]:
        """
        Returns the hit as a row in a table.
        :param report_section: Section is passed to save the alignments
        :param sub_dir: Specific subdirectory of the base directory to store report files
        :return: Table row
        """
        pass

    @abc.abstractmethod
    def get_table_column_names(self) -> List[str]:
        """
        Returns the table column names.
        :return: Table column names
        """
        pass

    @abc.abstractmethod
    def get_html_column_names(self) -> List[str]:
        """
        Returns the HTML column names.
        :return: HTML column names
        """
        pass

    @abc.abstractmethod
    def is_perfect_hit(self) -> bool:
        """
        Returns true if this is a perfect hit.
        :return: True if perfect
        """
        pass
