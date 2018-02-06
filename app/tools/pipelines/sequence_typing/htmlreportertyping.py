import os

from app.components.filesystemhelper import FileSystemHelper
from app.components.html.htmlreportsection import HtmlReportSection
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.io.tooliovalue import ToolIOValue
from app.tools.tool import Tool


class HtmlReporterTyping(Tool):
    """
    Tool that creates HTML reports for the sequence typing pipeline.

    Input:
        - HTML: Path to the HTML file to write the report
        - DIR: Directory to store files that are included in the HTML report
        - Informs 'Scheme': Information about the scheme
        - VAL_Hits: Hits detected for each locus
    Output:
        - HTML: Path to the generated report
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('HTML Reporter', '0.1', camel)
        self._report_section = None
        self._output_folder = None
        self._sub_folder = None

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__initialize_report()
        if 'ST' in self._input_informs:
            self.__add_sequence_type()
        self.__add_output_table()
        self._report_section.add_paragraph('Last updated: {}'.format(self._input_informs['Scheme']['last_updated']))
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section)]

    def __initialize_report(self):
        """
        Initializes the HTML report.
        :return: None
        """
        self._report_section = HtmlReportSection(self._input_informs['Scheme']['title'], 3)
        self._sub_folder = os.path.join(
            'sequence_typing', FileSystemHelper.make_valid(self._input_informs['Scheme']['name']))

    def __add_sequence_type(self):
        """
        Adds the sequence type to the report.
        :return: None
        """
        profile = self._input_informs['ST']['sequence_type']
        self._report_section.add_table([[value for _, value in profile.metadata]], [key for key, _ in profile.metadata],
                                       [('class', 'data')])

    def __add_output_table(self):
        """
        Adds the output table with the detected alleles.
        :return: None
        """
        header = self._tool_inputs['VAL_Hits'][0].value.get_html_column_names()
        hits = [t.value for t in self._tool_inputs['VAL_Hits']]
        table_data = [hit.to_html_row(self._report_section, self._sub_folder) for hit in hits]
        self._report_section.add_table(table_data, header, [('class', 'data')])

        tsv_file = self._tool_inputs['TSV'][0].path
        relative_path = os.path.join(self._sub_folder, os.path.basename(tsv_file))
        self._report_section.add_file(tsv_file, relative_path)
        self._report_section.add_link_to_file("Download (TSV)", relative_path)

    def _check_input(self):
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'Scheme' not in self._input_informs:
            raise InvalidInputSpecificationError("Scheme information is required")
        if 'VAL_Hits' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Typing hit input is required")
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Tabular input (TSV) is required")
        super()._check_input()
