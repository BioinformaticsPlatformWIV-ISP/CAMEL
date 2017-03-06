import os.path

from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.io.tooliofile import ToolIOFile
from app.tools.qiime.qiime import Qiime


class QiimeJoinPairedEnds(Qiime):
    """
    This script takes forward and reverse Illumina reads and joins them using the method chosen.
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super(QiimeJoinPairedEnds, self).__init__('qiime_join_paired_ends', '1.9.1', camel)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTQ_PE key is required
        - No additional keys allowed
        - Only two files allowed per key
        :return: None
        """
        if 'FASTQ_PE' not in self._tool_inputs:
            raise InvalidInputSpecificationError('Invalid input files (keys) given for join_paired_ends: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidInputSpecificationError('Too many input keys given for join_paired_ends: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs['FASTQ_PE']) != 2:
            raise InvalidInputSpecificationError('Invalid number (!= 2) of files in each key given for \
                                                  join_paired_ends: {!r}'.format(self._tool_inputs))

    def _set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs['FASTQ'] = [ToolIOFile(os.path.join(self._folder, 'fastqjoin.join.fastq'))]
        self._tool_outputs['FASTQ_Unjoined'] = [ToolIOFile(os.path.join(self._folder, 'fastqjoin.un1.fastq')),
                                                ToolIOFile(os.path.join(self._folder, 'fastqjoin.un2.fastq'))]

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        input_string = '-f {}'.format(self._tool_inputs['FASTQ_PE'][0])
        input_string += ' -r {}'.format(self._tool_inputs['FASTQ_PE'][1])
        input_string += ' -o {}'.format(self._folder)
        return input_string
