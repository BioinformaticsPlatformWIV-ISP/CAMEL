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

    def __init__(self, camel: Camel, input_type='illumina') -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        :param input_type: define the input data type of the analysis, defaults to illumina
        """
        super().__init__('ConFindr reporter', '0.1', camel)
        self._input_type = input_type

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
        section = HtmlReportSection('ConFindr', subtitle=self._input_informs['confindr']['_name'])

        # add message for ONT input
        if self._input_type in ('ont', 'hybrid'):
            section.add_warning_message('BE AWARE: running confindr with ONT data input is experimental! Results should be interpreted with the necessary caution!')
        # Main table
        cell = HtmlTableCell('No', 'green') if \
            self._input_informs['confindr']['NumContamSNVs'] < 20 else HtmlTableCell('Yes', 'red')
        section.add_header('Output', 3)
        if self._input_informs['confindr']['Genus'] == 'Error processing sample':
            section.add_alert('Error processing sample, species might be missing from rMLST database', 'warning')
        else:
            section.add_table([
                [cell,
                 f"<i>{self._input_informs['confindr']['Genus']}</i>",
                 self._input_informs['confindr']['NumContamSNVs'],
                 f"{self._input_informs['confindr']['BasesExamined']:,}"
                 ]],
                ['Contaminated', 'Genus', 'Contaminated SNPs', 'Bases examined'],
                table_attributes=[('class', 'data')]
            )

        # Last update
        try:
            update_date = date(*[int(x) for x in self._input_informs['confindr']['DatabaseDownloadDate'].split('-')])
        except ValueError:
            update_date = None
        update_date_str = update_date.strftime('%d-%m-%Y') if update_date is not None else 'NA'
        section.add_paragraph(f"Last database update: {update_date_str}")
        section.add_paragraph(
            '<b>Note:</b> Samples with multiple detected genera are always classified as contaminated.')

        # Set output
        self._tool_outputs['HTML'] = [ToolIOValue(section)]
