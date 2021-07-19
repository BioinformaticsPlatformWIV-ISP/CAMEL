import os
import re

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
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

    def _check_input(self) -> None:
        """
        Check and set the input specification.
        :return: None
        """
        super(Picard, self)._check_input()

        self._set_input()

    def _set_input(self) -> None:
        """
        Set the input specification
        :return: None
        """
        self._input_string = "I="
        if len(self._tool_inputs['BAM']) > 1:
            self._input_string += ' I='.join(f.path for f in self._tool_inputs['BAM'])
        else:
            self._input_string += self._tool_inputs['BAM'][0].path

        if 'FASTA_REF' in self._tool_inputs:
            self._input_string += f" R={self._tool_inputs['FASTA_REF'][0].path}"

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

    def _set_output(self) -> None:
        """
        Set the output specification
        :return: None
        """
        self._tool_outputs['BAM'] = [ToolIOFile(os.path.join(self._folder, self._parameters['output'].value))]
        self._tool_outputs['METRICS'] = [
            ToolIOFile(os.path.join(self._folder, self._parameters['matrics_output'].value))]
