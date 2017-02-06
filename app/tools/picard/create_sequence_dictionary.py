import os

from app.components.files.fileutils import FileUtils
from app.tools.picard.picard import Picard


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
        super(CreateSequenceDictionary, self).__init__('Picard CreateSequenceDictionary', '2.6.0', camel)

        self.function_name = 'CreateSequenceDictionary'
        self.required_inputs = ['FASTA_REF']
        self.supported_inputs = []
        self.specific_parameters = ['output_ext']

    def _set_output(self):
        """
        Set the output specification
        :return: None
        """
        # NOTE that the index file should be generated under the same directory as the FASTA_REF file
        index_file_name = os.path.splitext(self._tool_inputs['FASTA_REF'][0].path)[0] + \
            self._parameters['output_ext'].value
        self._output_string = ' O={}'.format(index_file_name)

        # NOTE if dictionary file exists, it should be removed otherwise Picard run will fail
        if os.path.isfile(index_file_name):
            FileUtils.silent_remove(index_file_name)
