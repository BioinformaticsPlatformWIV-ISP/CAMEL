from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.pipelines.shigella.subspeciesdetector import SubspeciesDetector
from camel.app.tools.tool import Tool


class SubspeciesReporter(Tool):
    """
    This tool is used to create a report for the Shigella subspecies determination.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance.
        :return: None
        """
        super().__init__('Shigella: subspecies reporter', '0.1', camel)
        self._section = HtmlReportSection('Subspecies identification')
        self._sub_folder = 'subspecies_identification'

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._section.add_paragraph('Detected subspecies: <b>{}</b>'.format(
            self._input_informs['subspecies']['detected_subspecies']))
        self.__add_profiles_table()
        self.__add_hits_table()
        self._tool_outputs['VAL_HTML_subspecies'] = [ToolIOValue(self._section)]

    def __add_profiles_table(self) -> None:
        """
        Adds a table with the profiles.
        :return: None
        """
        self._section.add_header('Profiles', 4)
        self._section.add_table(
            SubspeciesDetector.PROFILES['data'], SubspeciesDetector.PROFILES['header'], [('class', 'data')])

    def __add_hits_table(self) -> None:
        """
        Adds a table with the detected hits.
        :return: None
        """
        self._section.add_header('Hits', 4)
        table_data = [h.value.to_html_row(self._section, self._sub_folder, ) for h in
                      self._tool_inputs['VAL_hits'] if h.value.locus not in ('ipaH', 'speG')]
        self._section.add_table(table_data, self._input_informs['columns'], [('class', 'data')])
