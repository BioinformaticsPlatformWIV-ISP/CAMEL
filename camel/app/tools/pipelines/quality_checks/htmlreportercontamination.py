from datetime import date

import os
import re

from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class HtmlReporterContamination(Tool):
    """
    HTML reporter for the contamination check step.
    """

    ABBREVIATIONS = {
        'U': 'Unclassified',
        'D': 'Domain',
        'K': 'Kingdom',
        'P': 'Phylum',
        'C': 'Class',
        'O': 'Order',
        'F': 'Family',
        'G': 'Genus',
        'S': 'Species'
    }

    TITLE = 'Contamination check'

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('HTML Reporter: Contamination', '0.1', camel)
        self._subfolder = 'contamination_check'
        self._report_section = None

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self._report_section = HtmlReportSection(
            HtmlReporterContamination.TITLE, subtitle=self._input_informs['kraken']['_name'])
        self.__add_species_table()
        self.__add_detailed_table(self._tool_inputs['TSV'][0].path)
        self.__add_krona_report()
        self.__add_warnings()
        self.__add_last_update(self._tool_inputs['DB'][0].path)
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
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Kraken report input (TSV) is required")
        super()._check_input()

    def __add_last_update(self, db_path):
        """
        Adds the date of the last update.
        :param db_path: Database path
        :return: None
        """
        m = re.match('.*kraken/(\\d{4})(\\d{2})(\\d{2})/abfhpv.*$', os.path.realpath(db_path))
        if m:
            d = date(int(m.group(1)), int(m.group(2)), int(m.group(3)))
            value = d.strftime('%d-%m-%Y')
        else:
            value = 'NA'
        self._report_section.add_paragraph('Last updated: {}'.format(value))

    def __add_species_table(self):
        """
        Adds a table containing the detected species and corresponding percentages.
        :return: None
        """
        header = ['Species', 'Percentage']
        expected_name, expected_perc = self._input_informs['species']['expected']
        expected_perc = '{:.2f}'.format(float(expected_perc))
        table_data = [
            [HtmlElement('th', 'Expected', [('colspan', 2)])],
            [HtmlTableCell('<i>{}</i>'.format(expected_name), 'green'), HtmlTableCell(expected_perc, 'green')],
            [HtmlElement('th', 'Contaminants', [('colspan', 2)])]
        ]
        if (len(self._input_informs['species']['contaminants_fail']) +
                len(self._input_informs['species']['contaminants_warn']) == 0):
            table_data.append([HtmlTableCell('None found', attributes=[('colspan', 2)])])
        for species_name, percentage in self._input_informs['species']['contaminants_fail']:
            table_data.append([
                HtmlTableCell('<i>{}</i>'.format(species_name), 'red'),
                HtmlTableCell(percentage, 'red')])
        for species_name, percentage in self._input_informs['species']['contaminants_warn']:
            table_data.append([
                HtmlTableCell('<i>{}</i>'.format(species_name), 'orange'),
                HtmlTableCell(percentage, 'orange')])
        self._report_section.add_table(table_data, header, [('class', 'data')])

    def __add_detailed_table(self, kraken_report_path: str) -> None:
        """
        Adds a table with the detailed information from the tabular KRAKEN output.
        :return: None
        """
        table_data = []
        with open(kraken_report_path) as handle:
            for line in handle.readlines():
                parts = line.strip().split('\t')
                percentage = float(parts[0])
                if percentage > 1:
                    # Add genus / species in italics
                    if parts[3] in ('G', 'S'):
                        name = '<i>{}</i>'.format(parts[5].strip())
                    else:
                        name = parts[5]
                    table_data.append(['{:.2f}'.format(percentage), HtmlReporterContamination.ABBREVIATIONS.get(
                        parts[3], '-'), name])
        header = ['Percentage', 'Level', 'Name']
        self._report_section.add_table(table_data, header, [('class', 'data')])

    def __add_krona_report(self):
        """
        Adds a download link to the Krona report.
        :return: None
        """
        relative_path = os.path.join(self._subfolder, 'krona_report.html')
        self._report_section.add_file(self._tool_inputs['HTML_Krona'][0].path, relative_path)
        self._report_section.add_link_to_file('Krona Report', relative_path)

    def __add_warnings(self):
        """
        Adds warnings to the report when there is possible contamination.
        :return: None
        """
        if len(self._input_informs['species']['contaminants_warn']) > 0:
            self._report_section.add_warning_message("The sample is possibly contaminated with: {}".format(
                ", ".join('<i>{}</i>'.format(species_name) for species_name, _ in
                          self._input_informs['species']['contaminants_warn'])))
        if len(self._input_informs['species']['contaminants_fail']) > 0:
            self._report_section.add_error_message("There is strong evidence of a contamination with: {}".format(
                ", ".join('<i>{}</i>'.format(species_name) for species_name, _ in
                          self._input_informs['species']['contaminants_fail'])))

    @staticmethod
    def generate_empty_section():
        """
        Returns a report that is used when this analysis is disabled.
        :return: Report section
        """
        section = HtmlReportSection(HtmlReporterContamination.TITLE)
        section.add_paragraph('Analysis disabled')
        return section
