from typing import List, Tuple

from camel.app.camel import Camel
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class HtmlReporterQualityChecks(Tool):
    """
    Tool to create HTML report for the quality checks pipeline.
    """

    COLOR_CODES = {
        'OK': 'green',
        'Warning': 'orange',
        'Skipped': 'yellow',
        'Failed': 'red'
    }

    def __init__(self, camel: Camel) -> None:
        """
        Initialize this tool.
        :param camel: CAMEL instance
        :return: None
        """
        super().__init__('Quality checks reporter', '0.1', camel)
        self._sub_folder = 'quality_control'
        self._section = HtmlReportSection('Quality checks')

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self.__add_metric_overview()
        explanations = [
            (qc_check['full_name'], qc_check['explanation']) for qc_check in self._input_informs['qc_checks'] if
            (qc_check.get('explanation') is not None) and (qc_check.get('ori') != 'rev')]
        if len(explanations) > 0:
            self.__add_explanation_metrics(explanations)
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_metric_overview(self) -> None:
        """
        Adds the table with the metrics.
        :return: None
        """
        header = ['Metric', 'Status', 'Value', 'Warning threshold', 'Fail threshold']
        data = [[
            inf['full_name'] + (f" ({inf['ori']}.)" if inf.get('ori') is not None else ''),
            HtmlTableCell(inf['status'], color=HtmlReporterQualityChecks.COLOR_CODES[inf['status']]),
            inf['fmt_string_value'].format(inf['value']) if inf['value'] is not None else 'NA',
            inf['fmt_string_value'].format(inf['threshold_warn']),
            inf['fmt_string_value'].format(inf['threshold_fail'])
        ] for inf in self._input_informs['qc_checks']]
        self._section.add_header(f'Overview - {self._input_informs.get("read_type", "illumina").capitalize()} reads', 3)
        self._section.add_table(data, header, [('class', 'data')])

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        if 'qc_checks' not in self._input_informs:
            raise InvalidInputSpecificationError("No QC checks informs found")
        super(HtmlReporterQualityChecks, self)._check_input()

    def __add_explanation_metrics(self, explanations: List[Tuple[str, str]]) -> None:
        """
        Adds an explanation of the quality checks.
        :param explanations: List of explanations
        :return: None
        """
        self._section.add_header('Additional information', 3)
        self._section.add_labeled_list(explanations)
