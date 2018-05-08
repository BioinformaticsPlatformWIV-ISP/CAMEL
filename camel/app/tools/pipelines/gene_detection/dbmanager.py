import json

import os


from camel.app.components.genedetection.mapping import Mapping
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class DBManager(Tool):
    """
    Tool that manages gene databases.

    INPUT:
    - DIR: Directory containing an indexed FASTA file, a db_metadata.txt and a mapping.txt file.

    OUTPUT:
    - FASTA: FASTA file from the input directory

    INFORMS:
    - title: Database title
    - name: Database name
    - last_updated: Date of the last database update
    - clustering_cutoff: Clustering cutoff of the database
    - metadata: Metadata for all sequences in the FASTA file (parsed JSON from headers).
    - mapping: Mapping of converted sequence names to original headers
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
        self.__add_informs(input_folder)

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

    def __add_informs(self, input_folder):
        """
        Adds the informs.
        :param input_folder: Input database folder
        :return: None
        """
        if not os.path.isfile(os.path.join(input_folder, 'db_metadata.txt')):
            raise FileNotFoundError("No database metadata found in: {}".format(input_folder))
        with open(os.path.join(input_folder, 'db_metadata.txt')) as handle:
            metadata = json.load(handle)
        self._informs['title'] = metadata['title']
        self._informs['last_updated'] = metadata['last_updated']
        self._informs['clustering_cutoff'] = metadata['clustering_cutoff']
        self._informs['name'] = metadata['name']
        try:
            self._informs['mapping'] = Mapping.parse(os.path.join(input_folder, 'mapping.txt'))
        except FileNotFoundError:
            raise ToolExecutionError("No mapping found in {}".format(input_folder))

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
