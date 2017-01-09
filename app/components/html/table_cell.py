from app.components.html.element import HtmlElement


class HtmlTableCell(HtmlElement):
    """
    Class that represents a cell in a HTML table.
    """

    def __init__(self, text, attributes=None, link=None):
        """
        Initialize a new HTML element.
        :param text: Text displayed in th8e cell
        :param attributes: Attributes
        :param link: Hyperlink of the text
        """
        self._text_value = text
        self._attributes = attributes
        self._link = link
        super(HtmlTableCell, self).__init__()

    @property
    def text_value(self):
        """
        Returns the displayed text.
        :return: Text value
        """
        return self._text_value

    def _generate_html(self):
        """
        Generates the HTML representation of this element.
        :return: HTML representation
        """
        if self._link is not None:
            with self._get_tag('td', self._attributes):
                with self._get_tag('a', [('href', self._link)]):
                    self._text(self._text_value)
        else:
            with self._get_tag('td', self._attributes):
                self._text(self._text_value)
        return self._doc.getvalue()
