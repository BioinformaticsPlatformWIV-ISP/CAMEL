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

    def _set_input(self) -> None:
        """
        Set the input specification
        :return: None
        """
        # input can be one or multiple BAM or SAM files
        if len(self._tool_inputs[self._main_input]) > 1:
            self._input_string += ''.join(f"I={f.path} " for f in self._tool_inputs[self._main_input])
        else:
            self._input_string += f'I={self._tool_inputs[self._main_input][0].path} '

        self._tool_inputs.remove(self._main_input)

        # Run parent function for other inputs, e.g. FASTA_REF (optional)
        super(MarkDuplicates, self)._set_input()

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
            ToolIOFile(os.path.join(self._folder, self._parameters['metrics_output'].value))]
