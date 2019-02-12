import os

from camel.app.components.html.htmlelement import HtmlElement
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
        'Pass': 'green',
        'Warn': 'yellow',
        'Fail': 'red'
    }

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: CAMEL instance
        :return: None
        """
        super().__init__('Quality checks reporter', '0.1', camel)
        self._sub_folder = 'quality_control'
        self._section = HtmlReportSection('Additional quality checks')

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__add_metrics_overview()
        self.__add_fastqc_checks()
        self.__add_explanation_fastqc_checks()
        self.__add_warnings()
        self.__add_errors()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._section)]

    def _check_input(self) -> None:
        """
        Checks if the input is valid.
        :return: None
        """
        if 'st_genes' not in self._input_informs:
            raise InvalidInputSpecificationError("No sequence typing info found")
        if 'coverage' not in self._input_informs:
            raise InvalidInputSpecificationError("No coverage info found")
        if 'mapping' not in self._input_informs:
            raise InvalidInputSpecificationError("No mapping info found")
        if 'additional_checks' not in self._input_informs:
            raise InvalidInputSpecificationError("No additional quality check info found")
        if 'quality_criteria' not in self._input_informs:
            raise InvalidInputSpecificationError("No quality criteria info found")
        super(HtmlReporterQualityChecks, self)._check_input()

    def __add_metrics_overview(self) -> None:
        """
        Adds an overview with the QC metrics.
        :return: None
        """
        self._section.add_header('Metrics', 3)
        self.__add_qc_metric_with_plot(
            'Median coverage',
            '{:.0f}X'.format(self._input_informs['coverage']['median_depth']),
            'PNG_cov'
        )
        self.__add_qc_metric_with_plot(
            '{} genes found'.format(self._input_informs['st_genes']['title']),
            self.__get_st_stats(),
            'PNG_st'
        )
        self.__add_qc_metric_with_plot(
            'Reads mapping back to assembly',
            '{}%'.format(self._input_informs['mapping']['stats_map_rate']),
            'PNG_mapping'
        )

    def __add_qc_metric_with_plot(self, title: str, value: str, img_key: str) -> None:
        """
        Adds a QC metric to the output report, along with the corresponding plot.
        :param title: Title
        :param value: Value
        :param img_key: Image tool input key
        :return: None
        """
        self._section.add_paragraph(f'<b>{title}</b>: {value}')
        img_path = self._section.add_file(
            self._tool_inputs[img_key][0].path, os.path.join(self._sub_folder, f'plot_{img_key.split("_")[-1]}.png'))
        self._section.add_html_object(
            HtmlElement('img', attributes=[('src', img_path), ('width', 480), ('class', 'bordered')]))

    def __add_fastqc_checks(self) -> None:
        """
        Adds the table containing the additional QC checks.
        :return: None
        """
        self._section.add_header('FastQC additional checks', 3)
        informs = self._input_informs['additional_checks']
        table_data = []
        for test_name, test_results in sorted(informs['tests'].items()):
            table_data.append([test_name] + [self.__get_test_status_cell(result) for result in test_results])
        header = ['Test', 'Forward', 'Reverse'] if len(informs['samples']) == 2 else ['Test', 'Result']
        self._section.add_table(table_data, header, [('class', 'data')])

    @staticmethod
    def __get_test_status_cell(result: str) -> HtmlTableCell:
        """
        Returns a table cell for the given test result.
        :param result: Test result
        :return: HTML table cell
        """
        return HtmlTableCell(result, color=HtmlReporterQualityChecks.COLOR_CODES[result])

    def __add_explanation_fastqc_checks(self) -> None:
        """
        Adds an explanation of the quality checks.
        :return: None
        """
        test_explanations = [
            ['Average quality score test', 'checks whether the average read quality is above a threshold.'],
            ['GC content test', 'checks if the detected GC content is close enough to the expected GC content for this '
                                'organism.'],
            ['Maximal N-fraction test', 'checks whether the maximal N fraction at any read position is below a '
                                        'threshold.'],
            ['Mean Q-score drop test', 'checks whether the average position in the reads where the mean Q-score drops '
                                       'below <b>28</b> is above the warning / fail threshold.'],
            ['Per base sequence content test', 'checks whether the difference between A-T and C-G is below a threshold '
                                               'at every position. The beginning of the reads are skipped, as the peaks'
                                               ' there can be "normal".'],
            ['Sequence length distribution test', 'checks if the fraction of short sequences is below a threshold. The '
                                                  'warning and fail thresholds are determined dynamically based on the '
                                                  'mode length of the raw input reads (<b>{}</b>), the warning '
                                                  'threshold is set to 66.67% percent of this value (<b>{:.0f}</b>), '
                                                  'the fail threshold as 40.00% (<b>{:.0f}</b>).'
                .format(self._input_informs['additional_checks']['mode_read_length'],
                        self._input_informs['additional_checks']['length_warn'],
                        self._input_informs['additional_checks']['length_fail'])]
        ]
        self._section.add_labeled_list(test_explanations)

    def __get_st_stats(self) -> str:
        """
        Returns the sequence typing stats.
        :return: Formatted string with the sequence typing stats
        """
        st_stats = self._input_informs['st_genes']
        return '{0:}/{1:} ({2:.2f}%)'.format(
            st_stats['hits_found'], st_stats['nb_of_loci'],
            100 * float(st_stats['hits_found']) / st_stats['nb_of_loci'])

    def __add_warnings(self) -> None:
        """
        Adds warnings to the report.
        :return: None
        """
        for warning in self._input_informs['quality_criteria']['warnings']:
            self._section.add_warning_message(warning)

    def __add_errors(self) -> None:
        """
        Adds the error messages to the report.
        :return: None
        """
        for error in self._input_informs['quality_criteria']['fails']:
            self._section.add_error_message('{}, pipeline results can be inaccurate'.format(error))
