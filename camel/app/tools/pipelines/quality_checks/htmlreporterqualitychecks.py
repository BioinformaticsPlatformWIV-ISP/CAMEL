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
        self.__add_explanation_fastqc_metrics()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def __add_metric_overview(self) -> None:
        """
        Adds the table with the metrics.
        :return: None
        """
        header = ['Metric', 'Status', 'Value', 'Warning threshold', 'Fail threshold']
        data = [[
            inf['full_name'] + (f" ({inf['ori']}.)" if inf['ori'] is not None else ''),
            HtmlTableCell(inf['status'], color=HtmlReporterQualityChecks.COLOR_CODES[inf['status']]),
            inf['fmt_string_value'].format(inf['value']) if inf['value'] is not None else 'NA',
            inf['fmt_string_value'].format(inf['threshold_warn']),
            inf['fmt_string_value'].format(inf['threshold_fail'])
        ] for inf in self._input_informs['qc_checks']]
        self._section.add_header('Overview', 3)
        self._section.add_table(data, header, [('class', 'data')])

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        if 'qc_checks' not in self._input_informs:
            raise InvalidInputSpecificationError("No QC checks informs found")
        super(HtmlReporterQualityChecks, self)._check_input()

    def __add_explanation_fastqc_metrics(self) -> None:
        """
        Adds an explanation of the quality checks.
        :return: None
        """
        parser_informs = self._input_informs['fastqc_parser']
        test_explanations = [
            ['Average quality score test',
             'checks if the average read quality is above the given threshold.'],
            ['GC content test',
             'checks if the detected GC content is close enough to the expected GC content for this organism '
             '(<b>{:.2f}%</b>).'.format(float(self._parameters['gc_content_ref'].value))],
            ['Maximal N-fraction test',
             'checks if the maximal N fraction at any read position is below the given threshold.'],
            ['Mean Q-score drop test',
             'checks whether the average position in the reads where the mean Q-score drops below <b>{}</b> is above '
             'the given threshold.'.format(parser_informs['params']['qscore_drop_pos']['threshold'])],
            ['Per base sequence content test',
             'checks if the difference between A-T and C-G is below the given threshold at every position. '
             'The first {} and last {} bases of the reads are skipped, as the peaks there can be caused by the library'
             'kit or trimming artifacts.'.format(
                 parser_informs['params']['per_b_seq_content']['skipped_start'],
                 parser_informs['params']['per_b_seq_content']['skipped_end'])],
            ['Sequence length distribution test',
             'checks if the median read length of the trimmed reads is below a threshold compared to the mode length '
             'of the raw input reads (<b>{}</b>)'.format(parser_informs['stats']['mode_read_length_raw'])]
        ]
        self._section.add_header('Explanation: FastQC additional metrics', 3)
        self._section.add_labeled_list(test_explanations)
