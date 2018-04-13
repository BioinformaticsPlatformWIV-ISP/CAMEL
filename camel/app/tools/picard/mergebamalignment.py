import re

from camel.app.tools.picard.picard import Picard


class MergeBamAlignment(Picard):

    """
    Class for Picard MergeBamAlignment function
    """

    def __init__(self, camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super(MergeBamAlignment, self).__init__('Picard MergeBamAlignment', '2.8.3', camel)

        self._function_name = 'MergeBamAlignment'
        self._supported_inputs = []
        self._required_inputs = ['BAM_UNMAPPED', 'BAM_ALIGNED', 'FASTA_REF']

    def _set_input(self):
        """
        Set required inputs specification
        :return: None
        """
        super(MergeBamAlignment, self)._set_input()

        self._input_string += " UNMAPPED={} ALIGNED={}".format(
            self._tool_inputs['BAM_UNMAPPED'][0].path,
            self._tool_inputs['BAM_ALIGNED'][0].path
        )

    def _set_informs(self):
        """
        Analyse the result of picard run and update tool.informs
        :return: None
        """
        for l in self.stdout.splitlines():
            m = re.search('Finished reading (\d+) total records from alignment SAM/BAM.', l)
            if m:
                self.informs['reads_total'] = m.group(1)
            m = re.search('Wrote (\d+) alignment records and (\d+) unmapped reads.', l)
            if m:
                if m.group(1):
                    self.informs['aligned_reads'] = m.group(1)
                if m.group(2):
                    self.informs['unmapped_reads'] = m.group(2)
