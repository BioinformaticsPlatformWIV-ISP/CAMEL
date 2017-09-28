import json
import logging
import os

from app.components.filesystemhelper import FileSystemHelper
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.error.toolexecutionerror import ToolExecutionError
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class LocusManager(Tool):
    """
    Loads metadata and FASTA file from a sequence typing locus directory.

    Input:
        DIR: Directory that contains an indexed FASTA file and a locus_metadata.tsv file with the metadata for that
             locus.

    Output:
        FASTA: FASTA file from the input directory

    Informs:
        The metadata from the metadata file is parsed and added unaltered to the informs.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super().__init__('Typing: Locus Manager', '0.1', camel)

    def _check_input(self):
        """
        Checks if the specified input is correct.
        :return: None
        """
        if 'DIR' not in self._tool_inputs:
            raise InvalidInputSpecificationError("DIR input is required")
        super()._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        locus_folder = self._tool_inputs['DIR'][0].path
        fasta_file = FileSystemHelper.get_file_with_extension(locus_folder, '.fasta')
        self._tool_outputs['FASTA'] = [ToolIOFile(fasta_file)]
        self._informs = self.__get_metadata(os.path.join(locus_folder, 'locus_metadata.txt'))

    @staticmethod
    def __get_metadata(metadata_file):
        """
        Retrieves the metadata from a JSON file.
        :param metadata_file: Metadata file
        :return: Parsed metadata
        """
        try:
            with open(metadata_file) as handle:
                return json.load(handle)
        except FileNotFoundError:
            message = "Locus metadata not found in '{}'".format(os.path.dirname(metadata_file))
            logging.error(message)
            raise ToolExecutionError(message)
