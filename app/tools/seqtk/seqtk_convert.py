import os

from app.io.tooliofile import ToolIOFile
from app.tools.seqtk.seqtk import Seqtk


class SeqtkConvert(Seqtk):
    """
    Class that converts fastq/fasta file using seqtk
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

    def _check_input(self):
        """
        Check and set self._input_files based on self._tool_inputs, this function supports different type of input than superclass
        :return: None
        """
        super(Seqtk, self)._check_input()

        if 'FASTQ' not in self._tool_inputs:
            raise KeyError(
                'Seqtk function {!r} required FASTQ file is not specified in tool_inputs: {}!'.format(self._function_name, self._tool_inputs))

        elif len(self._tool_inputs['FASTQ']) != 1:
            raise ValueError(
                "Seqtk function {} supports only one input file of FASTQ type got {!r}.".format(
                    self._function_name, self._tool_inputs['FASTQ']))

        self._input_string = self._tool_inputs['FASTQ'][0].path

    def _set_output(self):
        """
        Set the output specification
        :return: None
        """
        self._output_string = os.path.join(self._folder, self._parameters['output_file'].value)
        self._tool_outputs['FASTA'] = [ToolIOFile(self._output_string)]
