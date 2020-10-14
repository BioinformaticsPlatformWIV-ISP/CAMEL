import re

from camel.app.camel import Camel
from camel.app.tools.picard.picard import Picard


class MarkDuplicates(Picard):

    """
    Class for Picard MarkDuplicates function
    """

    def __init__(self, camel: Camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard MarkDuplicates', '2.23.3', camel)

        self._function_name = 'MarkDuplicates'
        self._supported_inputs = ['BAM']

    def _set_informs(self) -> None:
        """
        Analyse the result of picard run and update tool.informs
        :return: None
        """
        for line in self.stdout.splitlines():
            m = re.search(r'Read (\d+) records. (\d+) pairs never matched', line)
            if m:
                self.informs['reads_total'] = m.group(1)
                self.informs['pairs_unmatched'] = m.group(2)

            m = re.search(r'Marking (\d+) records as duplicates', line)
            if m:
                self.informs['duplicates_count'] = m.group(1)

            m = re.search(r'Found (\d+) optical duplicate clusters', line)
            if m:
                self.informs['optical_duplicate_clusters_count'] = m.group(1)
