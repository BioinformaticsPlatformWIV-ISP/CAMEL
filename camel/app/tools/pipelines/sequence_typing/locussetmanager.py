import json
import os

from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.tools.tool import Tool


class LocusSetManager(Tool):
    """
    This tool extracts the metadata from a scheme / locus set directory.

    Input:
        DIR: Directory with several subdirectories for each locus and a scheme_metadata.txt file.

    Informs:
        The scheme metadata is added directly to the informs.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super().__init__('Typing: Locus Set Manager', '0.1', camel)

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
        locus_set_folder = self._tool_inputs['DIR'][0].path
        metadata_path = os.path.join(locus_set_folder, 'scheme_metadata.txt')
        if not os.path.isfile(metadata_path):
            raise ToolExecutionError("No scheme metadata found in '{}'".format(metadata_path))
        with open(metadata_path) as handle:
            self._informs = json.load(handle)
