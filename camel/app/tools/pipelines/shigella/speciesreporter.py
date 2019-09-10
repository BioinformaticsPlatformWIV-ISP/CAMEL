from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.pipelines.shigella.speciesdetector import SpeciesDetector
from camel.app.tools.tool import Tool


class SpeciesReporter(Tool):
    """
    This tool is used to create a report for the Shigella species determination.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance.
        :return: None
        """
        super().__init__('Shigella: species reporter', '0.1', camel)
        self._section = HtmlReportSection('Species identification')
        self._sub_folder = 'subspecies_identification'

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._section.add_paragraph('Detected species: <b>{}</b>'.format(
            self._input_informs['species']['detected_species']))
        self.__add_overview_table()
        self.__add_profiles_table()
        self.__add_hits_table()
        self._tool_outputs['VAL_HTML_species'] = [ToolIOValue(self._section)]

    def __add_overview_table(self) -> None:
        """
        Adds an table with an overview of the markers.
        :return: None
        """
        table_data = [
            ['<i>ipaH</i> present:', 'Yes' if self._input_informs['species']['ipaH_present'] else 'No'],
            ['<i>speG</i> present:', '{} (depth: {}X)'.format(
                'Yes' if self._input_informs['species']['speG_present'] is True else 'No',
                self._input_informs['species']['speG_depth'])],
            ['<i>speG</i> TG missing:', 'Yes' if self._input_informs['species']['speG_indel_present'] else 'No']
        ]
        self._section.add_table(table_data, table_attributes=[('class', 'information')])
        self._section.add_paragraph(f"speG is considered present if the coverage of the region is at least "
                                    f"{self._input_informs['species']['speG_depth_threshold']}X.")

    def __add_profiles_table(self) -> None:
        """
        Adds a table with the profiles.
        :return: None
        """
        self._section.add_header('Profiles', 4)
        self._section.add_table(
            SpeciesDetector.PROFILES['data'], SpeciesDetector.PROFILES['header'], [('class', 'data')])

    def __add_hits_table(self) -> None:
        """
        Adds a table with the detected hits.
        :return: None
        """
        self._section.add_header('Hits', 4)
        table_data = [h.value.to_html_row(self._section, self._sub_folder, ) for h in
                      self._tool_inputs['VAL_hits'] if h.value.locus in ('ipaH', 'speG')]
        self._section.add_table(table_data, self._input_informs['columns'], [('class', 'data')])
