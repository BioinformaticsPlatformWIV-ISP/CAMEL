from pathlib import Path

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class AmpliGoneFasta2BedReporter(Tool):
    """
    Reporting class for the AmpliGone FASTA to BED tool.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes the tool.
        :param camel: CAMEL instance.
        """
        super().__init__('AmpliGone fasta2bed reporter', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'BED' not in self._tool_inputs:
            raise InvalidInputSpecificationError('BED input is required')
        if 'ampligone' not in self._input_informs:
            raise InvalidInputSpecificationError("AmpliGone informs are required ('ampligone')")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Create report section
        section = HtmlReportSection('Primer localization', 3, subtitle=self._input_informs['ampligone']['_name'])

        # Information table
        section.add_table([
            ['Primers file:', self._input_informs['ampligone']['fasta_primers']],
            ['Error rate:', self._input_informs['ampligone']['primer_mismatch_rate']],
        ], None, [('class', 'information')])

        # Stats table
        section.add_table([[
            str(len(self._input_informs['ampligone']['primers_in'])),
            str(sum(nb > 0 for _, nb in self._input_informs['ampligone']['primers_out'].items())),
        ]], ['Nb. primers', 'Nb. mapped primers'], [('class', 'data')])

        # Output and info section
        relative_path = Path('preprocess', 'primer_locations.bed')
        section.add_link_to_file('Download primer locations (BED)', relative_path)
        section.add_file(self._tool_inputs['BED'][0].path, relative_path)

        # Set the tool outputs
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
