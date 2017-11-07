import logging

from app.components.html.htmlelement import HtmlElement


class HtmlReportSection(HtmlElement):
    """
    This class can be used to create a section in the HTML report.
    """

    def __init__(self, title, level=2):
        """
        Initializes a report section.
        :param title: Section title
        :param level: Header level
        """
        super(HtmlReportSection, self).__init__('div', attributes=[('class', 'report_section')])
        if title is not None:
            self.add_header(title, level)
        self._files = []

    @property
    def files(self):
        """
        Returns the files that were added to this report.
        :return: Files
        """
        return self._files

    def add_file(self, input_file, relative_path):
        """
        Adds the file to the report.
        :param input_file: Input file
        :param relative_path: Relative path
        :return: None
        """
        logging.info("Adding file to report section: {}".format(relative_path))
        self._files.append((input_file, relative_path,))
