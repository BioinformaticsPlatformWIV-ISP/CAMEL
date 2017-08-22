import os

from app.io.tooliofile import ToolIOFile
from app.tools.seqtk.seqtk import Seqtk


class SeqtkConvert(Seqtk):

    """
    Class that converts fastq into fasta file using seqtk
    """

    def __init__(self, camel):
        """
        Initialize seqtk
        :param camel: Camel instance
        :return: None
        """
        super(SeqtkConvert, self).__init__('Seqtk Convert', '1.2', camel)

        self._function_name = 'Convert'
        self._supported_inputs = ['FASTQ']
        self._specific_parameters = ['output_file']

    def _set_input_string(self):
        """
        Set the input specification
        :return: input_string containing input specification
        """
        return self._tool_inputs['FASTQ'][0].path

    def _set_output(self):
        """
        Set the output specification
        :return: None
        """
        self._output_string = os.path.join(self._folder, self._parameters['output_file'].value)
        self._tool_outputs['FASTA'] = [ToolIOFile(self._output_string)]
