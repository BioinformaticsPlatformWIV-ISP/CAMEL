import abc


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

    @staticmethod
    def compose_allele_link_url(locus, allele_id):
        """
        """
        if allele_id not in ('-', '?'):
            return "http://bigsdb.pasteur.fr/perl/bigsdb/bigsdb.pl?db=pubmlst_listeria_seqdef_public&page=alleleInfo&locus={}&allele_id={}".format(locus, allele_id)
        else:
            None
