from yattag import Doc
from yattag import indent

from app.components.html.tag import HtmlTag
from app.components.html.element import HtmlElement


class HtmlHelper(object):
    """
    Class that generates HTML reports.
    """

    def __init__(self, filename):
        """
        Open the given HTML file for writing.
        :param filename: Filename
        """
        self.filename = filename
        self.doc, self.tag, self.text = Doc().tagtext()

    def get_tag(self, tag, attributes=None):
        """
        Returns a tag with the given attributes.
        :param tag: Tag
        :param attributes: HTML attributes
        :return: Tag
        """
        if attributes:
            return self.tag(tag, *attributes)
        else:
            return self.tag(tag)

    def initialize(self, title, css_file=None):
        """
        Initializes a HTML report.
        :param title: Title of the report
        :param css_file: CSS style
        :return: None
        """
        self.doc.asis('<!DOCTYPE HTML>')
        with self.tag('head'):
            with self.tag('title'):
                self.text(title)
            self.doc.stag('meta', charset='UTF-8')
            if css_file:
                with open(css_file) as css, self.tag('style', type="text/css"):
                    self.doc.asis(css.read())

    def close(self):
        """
        Close the current report
        :return: None
        """
        with open(self.filename, 'a') as html_file:
            html_file.write(indent(self.doc.getvalue()))

    def add_line_break(self):
        """
        Adds a line break.
        :return: None
        """
        self.doc.stag('br')

    def add_horizontal_line(self):
        """
        Adds a horizontal line.
        :return: None
        """
        self.doc.stag('hr')

    def add_header(self, text, level, attributes=None):
        """
        Adds a header.
        :param text: Header text
        :param level: Header level (1-6)
        :param attributes: Additional HTML arguments (key-value)
        :return: None
        """
        with self.get_tag('h%s' % str(level), attributes):
            self.text(text)

    def add_paragraph(self, text, attributes=None):
        """
        Adds a text paragraph
        :param text: Paragraph text
        :param attributes: Additional HTML arguments (key-value)
        :return: None
        """
        with self.get_tag('p', attributes):
            self.text(text)

    def add_list(self, elements, ordered=False, attributes=None):
        """
        Adds a list.
        :param elements: List elements
        :param ordered: If True an ordered list is added
        :param attributes: Additional HTML arguments (key-value)
        :return: None
        """
        if ordered:
            list_tag = 'ol'
        else:
            list_tag = 'ul'

        with self.get_tag(list_tag, attributes):
            for value in elements:
                with self.tag('li'):
                    if isinstance(value, HtmlTag):
                        self.doc.asis(value.html_code)
                    else:
                        self.text(value)

    def add_labeled_list(self, rows, ordered=False, attributes=None):
        """
        Adds a labeled list.
        :param rows: List of label-text pairs
        :param ordered: If True an ordered list is added
        :param attributes: Additional HTML arguments (key-value)
        :return: None
        """
        if ordered:
            list_tag = 'ol'
        else:
            list_tag = 'ul'

        with self.get_tag(list_tag, attributes):
            for label, text in rows:
                with self.tag('li'):
                    with self.tag('b'):
                        self.text('{}: '.format(label))
                    self.text(text)

    def add_link_to_file(self, link_text, file_):
        """
        Adds a link to a file.
        :param link_text: Text of the link
        :param file_: Location of the file
        :return: None
        """
        if file_ is None:
            self.text(link_text)
            with self.tag('em', klass='not-available'):
                self.text('Not available')
        else:
            with self.tag('a', href=file_):
                self.text(link_text)
        self.add_line_break()

    def add_table(self, data, column_names=None, table_attributes=None):
        """
        Adds a table.
        :param data: Table data
        :param column_names: Column names
        :param table_attributes: Attributes of the table
        :return: None
        """
        with self.get_tag('table', table_attributes):
            if column_names:
                self._add_table_header(column_names)
            for row in data:
                self._add_table_row(row)
        self.add_line_break()

    def _add_table_header(self, column_names):
        """
        Adds a header to a table.
        :param column_names: Names of the columns
        :return: None
        """
        with self.get_tag('tr'):
            for name in column_names:
                with self.get_tag('th'):
                    self.text(name)

    def _add_table_row(self, data):
        """
        Adds a row to a table.
        :param data: Row data
        :return: None
        """
        # TODO: Check
        with self.get_tag('tr'):
            for value in data:
                if isinstance(value, HtmlElement):
                    self.doc.asis(value.html_code)
                else:
                    with self.get_tag('td'):
                        self.text(value)

    def add_error_message(self, message):
        """
        Adds an error message.
        :param message: Message text
        :return: None
        """
        with self.tag('div', klass='alert-box error'):
            with self.tag('span'):
                self.text('error: ')
            self.text(message)

    def add_warning_message(self, message):
        """
        Adds a warning message.
        :param message: Message text
        :return: None
        """
        with self.tag('div', klass='alert-box warning'):
            with self.tag('span'):
                self.text('warning: ')
            self.text(message)

    @DeprecationWarning
    def add_raw_html(self, html_code):
        """
        Adds a piece of HTML code to the document
        :param html_code: HTML code
        :return: None
        """
        self.doc.asis(html_code)
