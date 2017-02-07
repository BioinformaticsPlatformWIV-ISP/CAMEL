import re

from app.tools.picard.picard import Picard


class FastqToSam(Picard):
    """
    Class for Picard FastqToSam function
    """
    SAMPLE_NAME = 'sampleA'

    def __init__(self, camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super(FastqToSam, self).__init__('Picard FastqToSam', '2.6.0', camel)

        self._function_name = 'FastqToSam'

    def _check_input(self):
        """
        Check and set the input specification
        :return: None
        """
        if 'FASTQ_PE' in self._tool_inputs:
            self._input_string = 'FASTQ={} FASTQ2={}'.format(
                self._tool_inputs['FASTQ_PE'][0].path,
                self._tool_inputs['FASTQ_PE'][1].path
            )
        elif 'FASTQ_SE' in self._tool_inputs:
            self._input_string = 'FASTQ={}'.format(self._tool_inputs['FASTQ_SE'][0].path)
        else:
            raise KeyError('Picard FastqToSam requires FASTQ_SE or FASTQ_PE input.')

        if 'SAMPLE_NAME' in self._tool_inputs:
            self._input_string += " SM={0} RG={0}".format(self._tool_inputs['SAMPLE_NAME'][0].value)
        else:
            self._input_string += " SM={0} RG={0}".format(self.SAMPLE_NAME)

    def _set_inform(self):
        """
        Analyse the result of picard run and update tool.informs
        :return: None
        """
        for l in self.stdout.split('\n'):
            m = re.search('Auto-detected quality format as: ([a-zA-Z]+)', l)
            if m:
                if m.group(1) == 'Standard':
                    self.informs['quality_encoding'] = 'phred33'
                elif m.group(1) == 'Illumina':
                    self.informs['quality_encoding'] = 'phred64'
                elif m.group(1) == 'Solexa':
                    self.informs['quality_encoding'] = 'solexa66'
                else:
                    self.informs['quality_encoding'] = '[UNKNOWN]{}'.format(m.group(1))

            m = re.search('Processed (\d+) fastq reads', l)
            if m:
                self.informs['reads_processed'] = m.group(1)
