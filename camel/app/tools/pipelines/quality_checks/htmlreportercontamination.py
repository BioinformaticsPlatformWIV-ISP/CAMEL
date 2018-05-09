import os

from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class HtmlReporterContamination(Tool):

    """
    HTML reporter for the contamination check step.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('HTML Reporter: Contamination', '0.1', camel)
        self._subfolder = 'contamination_check'
        self._report_section = HtmlReportSection("Contamination Check")

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__add_species_table()
        self.__add_krona_report()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section)]

    def _check_input(self):
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'HTML_Krona' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Krona report input is required.")
        if 'species' not in self._input_informs:
            raise InvalidInputSpecificationError("Species input is required")
        super()._check_input()

    def __add_species_table(self):
        """
        Adds a table containing the detected species and corresponding percentages.
        :return: None
        """
        header = ['Species', 'Percentage']
        table_data = [[HtmlTableCell(v, 'green') for v in self._input_informs['species']['expected']]]
        for species_name, percentage in self._input_informs['species']['contaminants_fail']:
            table_data.append([HtmlTableCell(v, 'red') for v in (species_name, percentage)])
        for species_name, percentage in self._input_informs['species']['contaminants_warn']:
            table_data.append([HtmlTableCell(v, 'orange') for v in (species_name, percentage)])
        table_data.append(HtmlTableCell(v, 'grey') for v in self._input_informs['species']['unclassified'])
        self._report_section.add_table(table_data, header, [('class', 'data')])

    def __add_krona_report(self):
        """
        Adds a download link to the Krona report.
        :return: None
        """
        relative_path = os.path.join(self._subfolder, 'krona_report.html')
        self._report_section.add_file(self._tool_inputs['HTML_Krona'][0].path, relative_path)
        self._report_section.add_link_to_file('Krona Report', relative_path)
        # CURRENTLY DISABLED: test show 'krona_out.html.files' only contain reads to nodes association data. Not required for showing Krona figures.
        #
        # Add Krona HTML related files under /krona_out.html.files
        # krona_html_files_dir = os.path.join(
        #     os.path.dirname(self._tool_inputs['HTML_Krona'][0].path), 'krona_out.html.files')
        # for _file in os.listdir(krona_html_files_dir):
        #     file_relative_path = os.path.join(self._subfolder, 'krona_out.html.files', os.path.basename(_file))
        #     self._report_section.add_file(os.path.join(krona_html_files_dir, _file), file_relative_path)
