import os

from app.components.filesystemhelper import FileSystemHelper
from app.io.tooliodirectory import ToolIODirectory
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class GeneDBManager(Tool):
    """
    Tool that manages gene databases.

    INPUT:
    - DIR: Directory / symlink to a directory with the following structure:
        '../{DATABASE_NAME}/{DATE:YYYY/MM/DD}'

    OUTPUT:
    - FASTA: FASTA file found in the input directory

    INFORMS:
    - last_updated: Last update date from the database (derived from directory name)
    - database_name: Database name (derived from parent directory)
    """

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Resistance Characterization: Gene Database Manager', '0.1', camel)

    def _execute_tool(self):
        """
        Runs this tool.
        :return: None
        """
        input_folder = self._tool_inputs['DIR'][0].path
        self._tool_outputs['FASTA'] = [GeneDBManager.__get_database_file(input_folder)]
        self.__add_informs(input_folder)

    def _check_input(self):
        """
        Checks whether the input is correct.
        :return: None
        """
        if 'DIR' not in self._tool_inputs:
            raise ValueError("No 'DIR' input found.")
        if not isinstance(self._tool_inputs['DIR'][0], ToolIODirectory):
            raise ValueError("'{}' is not a directory".format(self._tool_inputs['DIR'][0]))
        super(GeneDBManager, self)._check_input()

    def __add_informs(self, input_folder):
        """
        Adds the informs.
        :return: None
        """
        name, last_updated = self.__get_database_info(input_folder)
        self._informs['name'] = name
        self._informs['last_updated'] = last_updated

    @staticmethod
    def __get_database_info(input_folder):
        """
        Returns the info on the locus.
        :param input_folder: Input folder
        :return: List of 'name' and 'last_updated' of the locus
        """
        dir_name = os.path.basename(os.path.realpath(input_folder))
        last_updated = '-'.join([dir_name[6:8], dir_name[4:6], dir_name[0:4]])
        name = os.path.split(os.path.dirname(input_folder))[-1]
        return name, last_updated

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
