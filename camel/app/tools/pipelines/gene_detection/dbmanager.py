import json
import os

from Bio import SeqIO

from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class DBManager(Tool):
    """
    Tool that manages gene databases.

    INPUT:
    - DIR: Directory containing an indexed FASTA file and a last_update file containing the date of the last database
        update.

    OUTPUT:
    - FASTA: FASTA file from the input directory

    INFORMS:
    - last_updated: Last update date from the database
    - name: Database name
    - metadata: Metadata for all sequences in the FASTA file (parsed JSON from headers).
    """

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Gene Detection: DB Manager', '0.1', camel)

    def _execute_tool(self):
        """
        Runs this tool.
        :return: None
        """
        input_folder = self._tool_inputs['DIR'][0].path
        fasta_file = DBManager.__get_database_file(input_folder)
        self._tool_outputs['FASTA'] = [fasta_file]
        self.__add_informs(input_folder, fasta_file)

    def _check_input(self):
        """
        Checks whether the input is correct.
        :return: None
        """
        if 'DIR' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No 'DIR' input found.")
        if not isinstance(self._tool_inputs['DIR'][0], ToolIODirectory):
            raise InvalidInputSpecificationError("'{}' is not a directory".format(self._tool_inputs['DIR'][0]))
        super(DBManager, self)._check_input()

    def __add_informs(self, input_folder, fasta_file):
        """
        Adds the informs.
        :param input_folder: Input database folder
        :param fasta_file: FASTA file
        :return: None
        """
        name, last_updated = self.__get_database_info(input_folder)
        self._informs['name'] = name
        self._informs['last_updated'] = last_updated
        self._informs['metadata'] = self.__parse_metadata(fasta_file)

    @staticmethod
    def __get_database_info(input_folder):
        """
        Returns the info on the locus.
        :param input_folder: Input folder
        :return: List of 'name' and 'last_updated' of the locus
        """
        try:
            with open(os.path.join(input_folder, 'last_update')) as handle:
                last_update = handle.read().strip()
        except FileNotFoundError:
            last_update = 'NA'
        name = os.path.basename(input_folder)
        return name, last_update

    @staticmethod
    def __parse_metadata(fasta_file):
        """
        Parses the FASTA file metadata.
        :param fasta_file: FASTA file
        :return: Metadata
        """
        metadata = {}
        with open(fasta_file.path) as handle:
            for seq in SeqIO.parse(handle, 'fasta'):
                metadata_str = ' '.join(seq.description.split(' ')[1:])
                metadata[seq.id] = json.loads(metadata_str)
        return metadata

    @staticmethod
    def __get_database_file(folder):
        """
        Returns the FASTA file from the given folder.
        :return: Database file
        """
        try:
            return ToolIOFile(FileSystemHelper.get_file_with_extension(folder, '.fasta'))
        except IOError:
            raise IOError("Cannot retrieve FASTA file from folder: '{}'".format(folder))
