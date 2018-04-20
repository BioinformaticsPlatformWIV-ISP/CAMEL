import logging
import unittest

import os
import tempfile

from camel.app.camel import Camel
from camel.app.components.html.htmlexpandablediv import HtmlExpandableDiv
from camel.app.components.html.htmlreport import HtmlReport
from camel.resources import CSS_STYLE
from camel.resources.javascript import JQUERY_SRC


class TestReadTrimming(unittest.TestCase):
    """
    Tests the HtmlReport and related classes.
    """

    camel = Camel()
    running_dir = None

    THREADS = 8

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(None, 'camel_', '/scratch/temp')

    def test_html_report(self):
        """
        Tests the creation and export of a simple HtmlReport.
        :return: None
        """
        report_path = os.path.join(self.running_dir, 'report.html')
        report = HtmlReport(report_path, self.running_dir)
        report.initialize("Test report", CSS_STYLE)
        report.add_header("Test report", 1)
        report.add_paragraph("This  is the report content")
        report.save()
        logging.info("Report saved in: {}".format(report_path))
        self.assertGreater(os.path.getsize(report_path), 0)

    def test_report_with_js(self):
        """
        Tests the creation of a report that includes Javascript.
        :return: None
        """
        report_path = os.path.join(self.running_dir, 'report.html')
        report = HtmlReport(report_path, self.running_dir, include_js=[JQUERY_SRC])
        report.initialize("Test report", CSS_STYLE)
        report.add_header("Test report", 1)
        div = HtmlExpandableDiv('large-table', 'Table')
        table_data = [['row', i] for i in range(0, 40)]
        div.add_table(table_data, ['Col 1', 'Col 2'], [('class', 'data')])
        report.add_html_object(div)
        report.save()
        logging.info("Report saved in: {}".format(report_path))
        self.assertGreater(os.path.getsize(report_path), 0)


if __name__ == '__main__':
    unittest.main()
