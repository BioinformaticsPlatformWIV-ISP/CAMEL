import re

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.tools.picard.picard import Picard


class FastqToSam(Picard):

    """
    Class for Picard FastqToSam function
    """

    def __init__(self, camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super(FastqToSam, self).__init__('Picard FastqToSam', '2.8.3', camel)

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
            InvalidInputSpecificationError('Picard FastqToSam requires FASTQ_SE or FASTQ_PE input.')

    def _set_input(self):
        """
        Set input for the picard function
        :return: None
        """
        if 'SAMPLE_NAME' in self._tool_inputs:
            # if SAMPLE_NAME specified, it will replace the default values of parameters: RG_sample_name, RG_name in DB
            self._specific_parameters = ['RG_sample_name', 'RG_name']
            self._input_string += " SM={0} RG={0}".format(self._tool_inputs['SAMPLE_NAME'][0].value)

    def _set_informs(self):
        """
        Analyse the result of picard run and update tool.informs
        :return: None
        """
        for l in self.stdout.splitlines():
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
