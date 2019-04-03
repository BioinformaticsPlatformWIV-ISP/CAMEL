import logging
import os

from camel.app.components.files.fileutils import FileUtils
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.picard.picard import Picard


class CreateSequenceDictionary(Picard):

    """
    Class for Picard CreateSequenceDictionary function
    """

    def __init__(self, camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Picard CreateSequenceDictionary', '2.8.3', camel)

        self._function_name = 'CreateSequenceDictionary'
        self._required_inputs = ['FASTA_REF']
        self._supported_inputs = []
        self._specific_parameters = ['output_ext', 'symlink']
        self._fasta_file = None

    def _set_output(self):
        """
        Set the output specification
        :return: None
        """
        # NOTE that the index file should be generated under the same directory as the FASTA_REF file
        index_file_name = os.path.splitext(self._fasta_file)[0] + self._parameters['output_ext'].value
        # NOTE if dictionary file exists, it should be removed otherwise Picard run will fail
        if os.path.isfile(index_file_name):
            FileUtils.silent_remove(index_file_name)
        self._output_string = ' O={}'.format(index_file_name)

        self._tool_outputs['FASTA_REF'] = [ToolIOFile(self._fasta_file)]

    def _set_input(self):
        """
        Function to set required and optional inputs in self._input_string
        :return: None
        """
        if 'symlink' in self._parameters:
            logging.info("Creating a symlink for the FASTA_REF input to generate SequenceDictionary locally.")
            self._fasta_file = self.__symlink_input()
        else:
            self._fasta_file = self._tool_inputs['FASTA_REF'][0].path

        logging.debug("symlink {}".format(self._fasta_file))

        self._input_string += " R={}".format(self._fasta_file)

    def __symlink_input(self):
        """
        Creates a symlink for the input.
        :return: Path to the symlink of the input
        """
        symlink_location = os.path.join(self._folder, self._tool_inputs['FASTA_REF'][0].basename)
        if not os.path.exists(symlink_location):
            os.symlink(self._tool_inputs['FASTA_REF'][0].path, symlink_location)

        return symlink_location
