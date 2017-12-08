import logging
import os

from app.components.filesystemhelper import FileSystemHelper
from app.components.html.htmlreportsection import HtmlReportSection
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.io.tooliovalue import ToolIOValue
from app.tools.tool import Tool


class HtmlReporterGeneDetection(Tool):

    """
    Tool that creates HTML reports for the resistance characterization pipeline.
    """

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Gene detection: reporter', '0.1', camel)
        self._subfolder = None
        self._report_section = None

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        db_name = self._input_informs['db_info']['name']
        self._report_section = HtmlReportSection(db_name, 3)
        self._subfolder = os.path.join('gene_detection', FileSystemHelper.make_valid(db_name))

        # Add output table
        if len(self._tool_inputs['VAL_Hits']) == 0:
            self._report_section.add_paragraph('No hits found.')
        else:
            header = self._tool_inputs['VAL_Hits'][0].value.get_html_column_names()
            table_data = [hit.to_html_row(self._report_section, self._subfolder) for hit in [
                t.value for t in self._tool_inputs['VAL_Hits']]]
            self._report_section.add_table(table_data, header, [('class', 'data')])

            # Add a download link to the TSV file
            tsv_file = self._tool_inputs['TSV'][0].path
            relative_path = os.path.join(self._subfolder, 'genes-{}-{}.tsv'.format(
                FileSystemHelper.make_valid(self._input_informs['db_info']['name']),
                FileSystemHelper.make_valid(self._tool_inputs['SAMPLE_NAME'][0].value)))
            self._report_section.add_file(tsv_file, relative_path)
            self._report_section.add_link_to_file("Download (TSV)", relative_path)

        # Add database information
        self._report_section.add_paragraph('Last updated: {}'.format(self._input_informs['db_info'].get(
            'last_updated', '{LAST_UPDATED}')))
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section)]

    def _check_input(self):
        """
        Checks if the input is valid.
        :return: None
        """
        if 'db_info' not in self._input_informs:
            raise InvalidInputSpecificationError("No database info found")
        if 'VAL_Hits' not in self._tool_inputs:
            logging.warning("No blast hits found")
        if 'SAMPLE_NAME' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Sample name input is required")
        super(HtmlReporterGeneDetection, self)._check_input()
