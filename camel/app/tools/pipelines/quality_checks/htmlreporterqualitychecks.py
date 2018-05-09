from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.components.html.htmltablecell import HtmlTableCell
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
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
        super().__init__('HTML Reporter', '0.1', camel)
        self.__sub_folder = 'quality_control'
        self._report_section = HtmlReportSection('Additional Quality Checks')

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        self.__add_additional_checks_section()
        self.__add_fastqc_checks()
        self.__add_explanation_fastqc_checks()
        self.__add_warnings()
        self.__add_errors()
        self._tool_outputs['VAL_HTML'] = [ToolIOValue(self._report_section)]

    def _check_input(self):
        """
        Checks if the input is valid.
        :return: None
        """
        if 'cgmlst' not in self._input_informs:
            raise InvalidInputSpecificationError("No cgMLST info found")
        if 'coverage' not in self._input_informs:
            raise InvalidInputSpecificationError("No coverage info found")
        if 'mapping' not in self._input_informs:
            raise InvalidInputSpecificationError("No mapping info found")
        if 'additional_checks' not in self._input_informs:
            raise InvalidInputSpecificationError("No additional quality check info found")
        if 'quality_criteria' not in self._input_informs:
            raise InvalidInputSpecificationError("No quality criteria info found")
        super(HtmlReporterQualityChecks, self)._check_input()

    def __add_additional_checks_section(self):
        """
        Adds the additional quality checks section.
        :return: None
        """
        data = [
            ['Median coverage:', self._input_informs['coverage']['median_depth']],
            ['cgMLST genes found:', self.__get_cgmlst_stats()],
            ['Reads mapping back to assembly:', '{}%'.format(self._input_informs['mapping']['stats_map_rate'])]
        ]
        self._report_section.add_table(data, table_attributes=[('class', 'information')])

    def __add_fastqc_checks(self):
        """
        Adds the table containing the additional QC checks.
        :return: None
        """
        informs = self._input_informs['additional_checks']
        table_data = []
        for test_name, test_results in sorted(informs['tests'].items()):
            table_data.append([test_name] + [self.__get_test_status_cell(result) for result in test_results])
        header = ['Test', 'Forward', 'Reverse']
        self._report_section.add_table(table_data, header, [('class', 'data')])

    @staticmethod
    def __get_test_status_cell(result):
        """
        Returns a table cell for the given test result.
        :param result: Test result
        :return: HTML table cell
        """
        return HtmlTableCell(result, color=HtmlReporterQualityChecks.COLOR_CODES[result])

    def __add_explanation_fastqc_checks(self):
        """
        Adds an explanation of the quality checks.
        :return: None
        """
        informs = self._input_informs['additional_checks']
        test_explanations = [
            ['Average quality score test', 'checks whether the average read quality is above a threshold.'],
            ['GC content test', 'checks if the detected GC content is close enough to the expected GC content for this organism.'],
            ['Maximal N-fraction test', 'checks whether the maximal N fraction at any read position is below a threshold.'],
            ['Mean Q-score drop test', 'checks whether the average position in the reads where the mean Q-score drops '
                                       'below <b>28</b> is above the warning / fail threshold.'],
            ['Per base sequence content test', 'checks whether the difference between A-T and C-G is below a threshold '
                                               'at every position. The beginning of the reads are skipped, as the peaks'
                                               ' there can be "normal".'],
            ['Sequence length distribution test', 'checks if the fraction of short sequences is below a threshold. The '
                                                  'warning and fail thresholds are determined dynamically based on the '
                                                  'length of the raw input reads (<b>{}</b>), the warning threshold is '
                                                  'set to 66.67% percent of this value (<b>{:.0f}</b>), the fail '
                                                  'threshold as 40.00% (<b>{:.0f}</b>).'.format(informs['max_read_length'], informs['length_warn'], informs['length_fail'])]
        ]
        self._report_section.add_labeled_list(test_explanations)

    def __get_cgmlst_stats(self):
        """
        Returns the cgMLST stats.
        :return: Formatted string with the cgMLST stats
        """
        cgmlst_stats = self._input_informs['cgmlst']
        return '{0:}/{1:} ({2:.2f}%)'.format(
            cgmlst_stats['hits_found'], cgmlst_stats['nb_of_loci'],
            100 * float(cgmlst_stats['hits_found']) / cgmlst_stats['nb_of_loci'])

    def __add_warnings(self):
        """
        Adds warnings to the report.
        :return: None
        """
        for warning in self._input_informs['quality_criteria']['warnings']:
            self._report_section.add_warning_message(warning)

    def __add_errors(self):
        """
        Adds the error messages to the report.
        :return: None
        """
        passed = True
        for error in self._input_informs['quality_criteria']['fails']:
            self._report_section.add_error_message('{}, pipeline aborted'.format(error))
            passed = False

        debug_mode = 'debug_mode' in self._parameters
        if not passed and not debug_mode:
            # TODO: Handle execution stops with Snakemake?
            pass
