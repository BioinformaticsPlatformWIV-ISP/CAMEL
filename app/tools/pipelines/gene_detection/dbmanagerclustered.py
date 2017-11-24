import json
import os

from app.components.filesystemhelper import FileSystemHelper
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.error.toolexecutionerror import ToolExecutionError
from app.io.tooliodirectory import ToolIODirectory
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class DBManagerClustered(Tool):
    """
    Tool that manages clustered gene databases.

    INPUT:
    - DIR: Directory containing an indexed FASTA file and a mapping file that contains the mapping to the original
        headers.

    OUTPUT:
    - FASTA: FASTA file from the input directory

    INFORMS:
    - Parsed mapping file
    """

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: Camel instance
        :return: None
        """
        super(DBManagerClustered, self).__init__('Gene Detection: DB Manager Clustered', '0.1', camel)

    def _execute_tool(self):
        """
        Runs this tool.
        :return: None
        """
        input_folder = self._tool_inputs['DIR'][0].path
        fasta_file = DBManagerClustered.__get_database_file(input_folder)
        self._tool_outputs['FASTA'] = [fasta_file]
        self._informs = self.__get_mapping(input_folder)

    def _check_input(self):
        """
        Checks whether the input is correct.
        :return: None
        """
        if 'DIR' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No 'DIR' input found.")
        if not isinstance(self._tool_inputs['DIR'][0], ToolIODirectory):
            raise InvalidInputSpecificationError("'{}' is not a directory".format(self._tool_inputs['DIR'][0]))
        super(DBManagerClustered, self)._check_input()

    @staticmethod
    def __get_mapping(input_folder):
        """
        Returns the mapping.
        :param input_folder: Input folder
        :return: Mapping as a dictionary
        """
        try:
            print('looking for mapping in {}'.format(input_folder))
            with open(os.path.join(input_folder, 'mapping.txt')) as handle:
                return json.load(handle)
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
