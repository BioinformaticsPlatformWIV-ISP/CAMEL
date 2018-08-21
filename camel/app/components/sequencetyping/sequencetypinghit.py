import abc
from typing import Optional


class SequenceTypingHit(metaclass=abc.ABCMeta):
    """
    This class represents a sequence typing hit.
    """

    def __init__(self, locus, allele_id):
        """
        Initializes the typing hit.
        :param locus: Locus
        :param allele_id: Allele id of the hit
        """
        self._locus = locus
        self._allele_id = allele_id
        self._allele_page_url_template = None
        # Listeria url template
        self.set_allele_page_url_template("http://bigsdb.pasteur.fr/perl/bigsdb/bigsdb.pl?db=pubmlst_listeria_seqdef_public&page=alleleInfo&locus={}&allele_id={}")

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
    def allele_id(self):
        """
        Returns the allele id.
        :return: Allele id
        """
        return self._allele_id

    @allele_id.setter
    def allele_id(self, allele_id):
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
        if self.allele_id.isdigit():
            return self._allele_page_url_template.format(self._locus, self.allele_id)
        else:
            return None

    @abc.abstractmethod
    def to_table_row(self):
        """
        Returns the hit as a row in a table.
        :return: Table row
        """
        pass

    @abc.abstractmethod
    def to_html_row(self, report_section, sub_dir=None):
        """
        Returns the hit as a row in a table.
        :param report_section: Section is passed to save the alignments
        :param sub_dir: Specific subdirectory of the base directory to store report files
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

    @abc.abstractmethod
    def is_perfect_hit(self):
        """
        Returns true if this is a perfect hit.
        :return: True if perfect
        """
        pass
