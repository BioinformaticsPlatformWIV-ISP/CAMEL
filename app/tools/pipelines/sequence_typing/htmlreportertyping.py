import os
import shutil

from app.components.filesystemhelper import FileSystemHelper
from app.components.html.htmlhelper import HtmlHelper
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class HtmlReporterTyping(Tool):
    """
    Tool that creates HTML reports for the sequence typing pipeline.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('Typing: HTML reporter', '0.1', camel)
        self._report = None
        self._output_folder = None
        self._sub_folder = None

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__initialize_report()
        self._report.add_header(self._input_informs['Scheme'].get('name', '{NAME}'), 3)
        if 'ST' in self._input_informs:
            self.__add_sequence_type()
        self.__add_output_table()
        self._report.add_paragraph('Last updated: {}'.format(self._input_informs['Scheme']['last_updated']))
        self._report.add_horizontal_line()
        self._report.close()
        self._tool_outputs['HTML'] = [ToolIOFile(self._report.filename)]

    def __initialize_report(self):
        """
        Initializes the HTML report.
        :return: None
        """
        self._report = HtmlHelper(self._tool_inputs['HTML'][0].path)
        self._output_folder = self._tool_inputs['DIR'][0].path
        self._sub_folder = os.path.join(
            'sequence_typing', FileSystemHelper.make_valid(self._input_informs['Scheme']['name']))
        if not os.path.isdir(os.path.join(self._output_folder, self._sub_folder, 'alignments')):
            os.makedirs(os.path.join(self._output_folder, self._sub_folder, 'alignments'))

    def __add_sequence_type(self):
        """
        Adds the sequence type to the report.
        :return: None
        """
        self._report.add_table([
            self._input_informs['ST']['metadata'].values()],
            self._input_informs['ST']['metadata'].keys(), [('class', 'data')])

    def __add_output_table(self):
        """
        Adds the output table with the detected alleles.
        :return: None
        """
        header = self._tool_inputs['VAL_Hits'][0].value.HTML_COLUMNS
        table_data = [hit.to_html_row(self._output_folder, self._sub_folder) for hit in [
            t.value for t in self._tool_inputs['VAL_Hits']]]
        self._report.add_table(table_data, header, [('class', 'data')])

        tsv_file = self._tool_inputs['TSV'][0].path
        shutil.copyfile(tsv_file, os.path.join(self._output_folder, self._sub_folder, os.path.basename(tsv_file)))
        self._report.add_link_to_file("Download (TSV)", os.path.join(self._sub_folder, os.path.basename(tsv_file)))

    def _check_input(self):
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'HTML' not in self._tool_inputs:
            raise InvalidInputSpecificationError("HTML input is required")
        if 'DIR' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Output directory is required")
        if 'Scheme' not in self._input_informs:
            raise InvalidInputSpecificationError("Scheme information is required")
        if 'VAL_Hits' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Typing hit input is required")
        super()._check_input()
