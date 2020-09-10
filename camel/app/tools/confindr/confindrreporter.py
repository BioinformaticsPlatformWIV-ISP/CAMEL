from datetime import date

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class ConFindrReporter(Tool):
    """
    Tool to generate reports for the ConFindr tool.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('ConFindr reporter', '0.1', camel)

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'confindr' not in self._input_informs:
            raise InvalidInputSpecificationError('ConFindr informs are required')
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        section = HtmlReportSection('ConFindr')

        # Main table
        cell = HtmlTableCell('No', 'green') if \
            self._input_informs['confindr']['ContamStatus'] is False else HtmlTableCell('Yes', 'red')
        section.add_header('Output', 3)
        section.add_table([
            [cell,
             f"<i>{self._input_informs['confindr']['Genus']}</i>",
             self._input_informs['confindr']['NumContamSNVs'],
             self._input_informs['confindr']['BasesExamined'],
             self._input_informs['confindr']['PercentContam'],
             self._input_informs['confindr']['PercentContamStandardDeviation']
             ]],
            ['Contaminated', 'Genus', 'Contaminated SNPs', 'Bases examined', 'Percent contaminated (%)',
             'Std. deviation (%)'],
            table_attributes=[('class', 'data')]
        )

        # Last update
        update_date = date(*[int(x) for x in self._input_informs['confindr']['DatabaseDownloadDate'].split('-')])
        section.add_paragraph(f"Last database update: {update_date.strftime('%d-%m-%Y')}")

        # Parameters
        section.add_header('Parameters', 3)
        section.add_table([
            ['Min. nb. of contaminant SNPs: ', 5],
            ['Min. nb. of bases to support multiple allele: ', 2],
        ], table_attributes=[('class', 'information')])
        section.add_paragraph(
            '<b>Note:</b> Samples with multiple detected genera are always classified as contaminated.')

        # Set output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
