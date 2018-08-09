import json
from typing import Tuple

import os

from camel.app.components.genedetection.mapping import Mapping
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
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
    - FASTA_clustered: Clustered FASTA file

    INFORMS:
    - title: Database title
    - name: Database name
    - last_updated: Date of the last database update
    - clustering_cutoff: Clustering cutoff of the database
    - mapping: Mapping of converted sequence names to original headers
    """

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Gene Detection: DB Manager', '0.1', camel)

    def _execute_tool(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        input_folder = self._tool_inputs['DIR'][0].path
        self.__set_database_files(input_folder)
        self.__add_informs(input_folder)

    def _check_input(self) -> None:
        """
        Checks whether the input is correct.
        :return: None
        """
        if 'DIR' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No 'DIR' input found.")
        if not isinstance(self._tool_inputs['DIR'][0], ToolIODirectory):
            raise InvalidInputSpecificationError("'{}' is not a directory".format(self._tool_inputs['DIR'][0]))
        super(DBManager, self)._check_input()

    def __add_informs(self, input_folder: str) -> None:
        """
        Adds the informs.
        :param input_folder: Input database folder
        :return: None
        """
        if not os.path.isfile(os.path.join(input_folder, 'db_metadata.txt')):
            raise FileNotFoundError("No database metadata found in: {}".format(input_folder))
        with open(os.path.join(input_folder, 'db_metadata.txt')) as handle:
            metadata = json.load(handle)
        self._informs.update(metadata)
        try:
            self._informs['mapping'] = self.__get_mapping(input_folder)
        except ToolExecutionError:
            raise FileNotFoundError(f'No mapping found in: {input_folder}')

    @staticmethod
    def __get_mapping(input_folder: str) -> Mapping:
        """
        Returns the mapping.
        :param input_folder: Input folder
        :return: Mapping as a dictionary
        """
        try:
            return Mapping.parse(os.path.join(input_folder, 'mapping.txt'))
        except FileNotFoundError:
            raise ToolExecutionError("No mapping found in {}".format(input_folder))

    @staticmethod
    def __get_database_info(input_folder: str) -> Tuple[str, str]:
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

    def __set_database_files(self, folder: str) -> None:
        """
        Returns the FASTA file from the given folder.
        :return: Database file
        """
        for f in os.listdir(folder):
            if f.endswith('.fasta') and 'clustered' in f:
                self._tool_outputs['FASTA_clustered'] = [ToolIOFile(os.path.join(folder, f))]
            elif f.endswith('.fasta'):
                self._tool_outputs['FASTA'] = [ToolIOFile(os.path.join(folder, f))]
