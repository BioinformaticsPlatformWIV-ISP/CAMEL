import abc
import os
import shutil
import tempfile

from camel.app.components.html.htmlreport import HtmlReport
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.tools.tool import Tool


class HtmlReporter(Tool):
    """
    Base class for creating HTML reports.
    """

    def __init__(self, camel, tool_name=None, tool_version=None):
        """
        Initializes this tool.
        :param camel: Camel instance
        :param tool_name: If specified this tool name is used, can be used to add parameters.
        :param tool_version: If specified this tool version is used, can be used to add parameters.
        """
        if tool_name is None:
            tool_name = 'HTML Reporter'
        if tool_version is None:
            tool_version = '0.1'
        super(HtmlReporter, self).__init__(tool_name, tool_version, camel)
        self._report = None
        self._output_folder = None

    @abc.abstractmethod
    def _create_report(self):
        """
        Creates the HTML report.
        :return: None
        """
        pass

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        if 'disable_html_output' not in self._parameters:
            self._report = HtmlReport(self._tool_inputs['HTML'][0].path)
        else:
            self._report = HtmlReport(tempfile.NamedTemporaryFile().name)
        self._output_folder = self._tool_inputs['DIR_HTML'][0].path
        self._create_report()
        self._report.save()
        self._tool_outputs['HTML'] = self._tool_inputs['HTML'].copy()

    def _check_input(self):
        """
        Checks if the required inputs were specified.
        :return: None
        """
        if 'HTML' not in self._tool_inputs:
            raise ValueError("No HTML input found")
        if 'DIR_HTML' not in self._tool_inputs:
            raise ValueError("No output directory found")
        if not isinstance(self._tool_inputs['DIR_HTML'][0], ToolIODirectory):
            raise ValueError("'DIR' input is not a directory")
        super(HtmlReporter, self)._check_input()

    def _save_file(self, file_, location):
        """
        Saves the file in the specified location (relative to the output directory).
        :param file_: File
        :param location: Relative location
        :return: None
        """
        if file_ is None:
            return None
        path_to_dir = os.path.join(self._output_folder, os.path.dirname(location))
        if not os.path.isdir(path_to_dir):
            os.makedirs(path_to_dir)
        shutil.copy(file_, os.path.join(path_to_dir, os.path.basename(location)))
        return location
