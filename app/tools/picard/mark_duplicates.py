import re

from app.tools.picard.picard import Picard


class MarkDuplicates(Picard):
    """
    Class for Picard MarkDuplicates function
    """

    def __init__(self, camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super(MarkDuplicates, self).__init__('Picard MarkDuplicates', '2.6.0', camel)

        self.function_name = 'MarkDuplicates'
        self.supported_inputs = ['BAM']

    def _set_inform(self):
        """
        Analyse the result of picard run and update tool.informs
        :return: None
        """
        for l in self.stdout.split('\n'):
            m = re.search('Read (\d+) records. (\d+) pairs never matched', l)
            if m:
                self.informs['reads_total'] = m.group(1)
                self.informs['pairs_unmatched'] = m.group(2)

            m = re.search('Marking (\d+) records as duplicates', l)
            if m:
                self.informs['duplicates_count'] = m.group(1)

            m = re.search('Found (\d+) optical duplicate clusters', l)
            if m:
                self.informs['optical_duplicate_clusters_count'] = m.group(1)
