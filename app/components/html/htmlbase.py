import abc
import bs4
import yattag
from bs4 import BeautifulSoup
from yattag import SimpleDoc


class HtmlBase(object):
    """
    Base class for HTML objects, contains common methods for classes that implement HTML objects.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        """
        Initializes an HTML base object.
        """
        self._doc = yattag.SimpleDoc()

    def get_tag(self, tag, attributes=None):
        """
        Returns a tag with the given attributes.
        :param tag: Tag
        :param attributes: HTML attributes
        :return: Tag
        """
        if attributes:
            return self._doc.tag(tag, *attributes)
        else:
            return self._doc.tag(tag)

    def add_header(self, text, level, attributes=None):
        """
        Adds a header.
        :param text: Header text
        :param level: Header level (1-6)
        :param attributes: Additional HTML arguments (key-value)
        :return: None
        """
        with self.get_tag('h{}'.format(str(level)), attributes):
            self.add_text(text)

    def add_paragraph(self, text, attributes=None):
        """
        Adds a text paragraph
        :param text: Paragraph text
        :param attributes: Additional HTML arguments (key-value)
        :return: None
        """
        with self.get_tag('p', attributes):
            self.add_text(text)

    def add_line_break(self):
        """
        Adds a line break.
        :return: None
        """
        self._doc.stag('br')

    def to_html(self):
        """
        Converts the HTML object to HTML code.
        :return: HTML code
        """
        return self._doc.getvalue()

    def add_table(self, data, column_names=None, table_attributes=None):
        """
        Adds a table.
        :param data: Table data
        :param column_names: Column names
        :param table_attributes: Attributes of the table
        :return: None
        """
        with self.get_tag('table', table_attributes):
            if column_names is not None:
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
                    self.add_text(name)

    def _add_table_row(self, data):
        """
        Adds a row to a table.
        :param data: Row data
        :return: None
        """
        with self.get_tag('tr'):
            for value in data:
                if isinstance(value, HtmlBase):
                    self._doc.asis(value.to_html())
                else:
                    with self.get_tag('td'):
                        self.add_text(value)

    def add_link_to_file(self, link_text, file_):
        """
        Adds a link to a file.
        :param link_text: Text of the link
        :param file_: Location of the file
        :return: None
        """
        if file_ is None:
            self._doc.text(link_text)
            with self._doc.tag('em', klass='not-available'):
                self._doc.text('Not available')
        else:
            with self._doc.tag('a', href=file_):
                self._doc.text(link_text)
        self.add_line_break()

    def add_error_message(self, message):
        """
        Adds an error message.
        :param message: Message text
        :return: None
        """
        with self._doc.tag('div', klass='alert-box error'):
            with self._doc.tag('span'):
                self._doc.text('error: ')
            self.add_text(message)

    def add_warning_message(self, message):
        """
        Adds a warning message.
        :param message: Message text
        :return: None
        """
        with self._doc.tag('div', klass='alert-box warning'):
            with self._doc.tag('span'):
                self._doc.text('warning: ')
            self.add_text(message)

    def add_labeled_list(self, rows, ordered=False, attributes=None):
        """
        Adds a labeled list.
        :param rows: List of label-text pairs
        :param ordered: If True an ordered list is added
        :param attributes: Additional HTML arguments (key-value)
        :return: None
        """
        list_tag = 'ol' if ordered else 'ul'
        with self.get_tag(list_tag, attributes):
            for label, text in rows:
                with self._doc.tag('li'):
                    with self._doc.tag('b'):
                        self._doc.text('{}: '.format(label))
                    self.add_text(text)

    def add_html_object(self, input_object):
        """
        Adds an HTML object.
        :param input_object: Input object
        :return: None
        """
        if not isinstance(input_object, HtmlBase):
            raise ValueError("{} is not an HTML object".format(input_object))
        self._doc.asis(input_object.to_html())

    def add_raw(self, html_code):
        """
        Adds raw HTML code, usage of this method is generally discouraged. Only when there is no other possibility.
        :param html_code: HTML code
        :return: None
        """
        self._doc.asis(html_code)

    def add_text(self, text):
        """
        Adds text to the this HTML tag.
        :return: None
        """
        if text is None:
            raise ValueError("Cannot add None as text to the HTML object")
        html = BeautifulSoup(str(text), 'html.parser')
        for part in html.contents:
            if isinstance(part, bs4.element.Tag):
                self._doc.asis(f'<{part.name}>')
                self._doc.text(part.text)
                self._doc.asis(f'</{part.name}>')
            else:
                self._doc.text(part)

    def add_module_header(self, text, id_=None):
        """
        Adds a sub header.
        :param text: Header text
        :param id_: Html id
        :return: None
        """
        attributes = [('class', 'sub_header')]
        if id_ is not None:
            attributes.append(('id', id_))
        with self.get_tag('div', attributes):
            with self.get_tag('h2'):
                self._doc.text(text)
