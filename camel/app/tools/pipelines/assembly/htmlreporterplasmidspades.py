from pathlib import Path

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class HTMLReporterPlasmidSpades(Tool):
    """
    Creates the HTML output for PlasmidSPAdes.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('HTML reporter: plasmidSPAdes', '1.0', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('FASTA input is required')
        if 'quast' not in self._input_informs:
            raise InvalidInputSpecificationError('QUAST informs are required')
        if 'spades' not in self._input_informs:
            raise InvalidInputSpecificationError('SPAdes informs are required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        section = HtmlReportSection('Assembly - plasmidSPAdes', 3, self._input_informs['spades']['_name'])
        quast_informs = self._input_informs['quast']
        section.add_table([
            ['N50:', '{:,}'.format(int(quast_informs['contig']['N50']))],
            ['Number of contigs:', '{:,}'.format(int(quast_informs['contig']['# contigs (>= 1000 bp)']))],
            ['Total length:', '{:,}'.format(int(quast_informs['genome']['Total length']))]
        ], table_attributes=[('class', 'information')])

        relative_path = Path('plasmidspades') / Path(self._tool_inputs['FASTA'][0].path).name
        section.add_link_to_file('Assembly - plasmidSPAdes (FASTA)', str(relative_path))
        section.add_file(self._tool_inputs['FASTA'][0].path, relative_path)

        self._tool_outputs['HTML'] = [ToolIOValue(section)]
