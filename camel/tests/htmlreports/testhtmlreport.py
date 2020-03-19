import logging
import unittest
from pathlib import Path

from camel.app.components.html.htmlexpandablediv import HtmlExpandableDiv
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.resources import CSS_STYLE
from camel.resources.javascript import JQUERY_SRC


class TestHtmlReporter(CamelTestSuite):
    """
    Tests the HtmlReport and related classes.
    """

    def test_html_report(self) -> None:
        """
        Tests the creation and export of a simple HtmlReport.
        :return: None
        """

        report_path = self.running_dir / 'report.html'
        report = HtmlReport(str(report_path), str(self.running_dir))
        report.initialize("Test report", CSS_STYLE)
        report.add_header("Test report", 1)
        report.add_paragraph("This  is the report content")
        report.save()
        logging.info("Report saved in: {}".format(report_path))
        self.assertGreater(report_path.stat().st_size, 0)

    def test_report_with_js(self) -> None:
        """
        Tests the creation of a report that includes Javascript.
        :return: None
        """
        report_path = self.running_dir / 'report.html'
        report = HtmlReport(str(report_path), str(self.running_dir), include_js=[JQUERY_SRC])
        report.initialize("Test report", CSS_STYLE)
        report.add_header("Test report", 1)
        div = HtmlExpandableDiv('large-table', 'Table')
        table_data = [['row', i] for i in range(0, 40)]
        div.add_table(table_data, ['Col 1', 'Col 2'], [('class', 'data')])
        report.add_html_object(div)
        report.save()
        logging.info("Report saved in: {}".format(report_path))
        self.assertGreater(report_path.stat().st_size, 0)

    def test_pipeline_header(self):
        """
        Tests adding a pipeline header to a report.
        :return: None
        """
        report_path = self.running_dir / 'report.html'
        report = HtmlReport(str(report_path), str(self.running_dir), include_js=[JQUERY_SRC])
        report.initialize("Test report", CSS_STYLE)
        report.add_pipeline_header('My <i>pipeline</i>')
        report.save()
        logging.info("Report saved in: {}".format(report_path))
        self.assertGreater(report_path.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
