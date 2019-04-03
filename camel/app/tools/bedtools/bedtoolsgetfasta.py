import os


from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.bedtools.bedtools import Bedtools


class BedtoolsGetFasta(Bedtools):

    """
    Bedtools GetFasta func class
    """
    DEFAULT_OUTPUT_NAME = 'bedtools_getfata_extracted_sequences.fa'

    def __init__(self, camel, tool_name='bedtools getfasta', version='2.25.0'):
        """
        Initialize a samtools tool.
        :param tool_name: Tool name
        :param version: Tool version
        :param camel: Camel instance
        :return: None
        """
        super().__init__(tool_name, version, camel)
        self._required_inputs = ['BED', 'FASTA']

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__set_output()
        self.__build_command()
        self._execute_command()

    def __build_command(self):
        """
        Builds the command.
        :return: None
        """
        self._command.command = ' '.join([
            self._tool_command,
            ' '.join(self._build_options()),
            '-bed {}'.format(self._tool_inputs['BED'][0].path),
            '-fi {}'.format(self._tool_inputs['FASTA'][0].path),
            '-fo {}'.format(self.DEFAULT_OUTPUT_NAME)
        ])

    def __set_output(self):
        """
        Sets the output of this tool.
        :return: None
        """
        self._tool_outputs['FASTA'] = [ToolIOFile(os.path.join(self._folder, self.DEFAULT_OUTPUT_NAME))]

    def _check_input(self):
        """
        Checks the input.
        :return: None
        """
        self._check_required_inputs()

        if len(self._tool_inputs['BED']) != 1:
            raise InvalidInputSpecificationError("Exactly one BED input file expected.")
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidInputSpecificationError("Exactly one FASTA input file expected.")

        super(BedtoolsGetFasta, self)._check_input()
