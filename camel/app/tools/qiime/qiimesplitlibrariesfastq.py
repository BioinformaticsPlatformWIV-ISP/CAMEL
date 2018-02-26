import os.path

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.qiime.qiime import Qiime


class QiimeSplitLibrariesFastq(Qiime):
    """
    Split_libraries_fastq demultiplexes fastq files if needed and performs several quality controls
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super(QiimeSplitLibrariesFastq, self).__init__('qiime_split_libraries_fastq', '1.9.1', camel)

    def _check_input(self):
        """
        Checks whether the given inputs are valid:
        - FASTA or FASTQ key is required
        - Only TSV_Map allowed as additional key
        - Only one file allowed per key
        :return: None
        """
        if 'FASTA' in self._tool_inputs == 'FASTQ' in self._tool_inputs:
            raise InvalidInputSpecificationError('Invalid input files (keys) given for '
                                                 'split_libraries_fastq: {!r}'.format(self._tool_inputs))
        for key, input_files in self._tool_inputs.items():
            if key not in ['FASTA', 'FASTQ', 'TSV_Map']:
                raise InvalidInputSpecificationError('Invalid input key given for split_libaries_fastq: {!r}'.format(self._tool_inputs))
            if len(input_files) != 1:
                raise InvalidInputSpecificationError('Invalid number (max = 1) of files in each key given for \
                                                     split_libraries_fastq: {!r}'.format(self._tool_inputs))

    def _set_output(self):
        """
        Sets the name of the output files
        :return: None
        """
        self._tool_outputs['TSV_Histogram'] = [ToolIOFile(os.path.join(self._folder, 'histograms.txt'))]
        self._tool_outputs['FASTA'] = [ToolIOFile(os.path.join(self._folder, 'seqs.fna'))]
        self._tool_outputs['LOG'] = [ToolIOFile(os.path.join(self._folder, 'split_library_log.txt'))]

    def _build_input_string(self):
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        if 'FASTA' in self._tool_inputs:
            input_string = '-i {}'.format(self._tool_inputs['FASTA'][0])
        else:
            input_string = '-i {}'.format(self._tool_inputs['FASTQ'][0])
        input_string += ' -o {}'.format(self._folder)
        return input_string
