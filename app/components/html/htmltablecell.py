from app.components.html.htmlelement import HtmlElement


class HtmlTableCell(HtmlElement):
    """
    This class represents a cell in a HTML table.
    """

    def __init__(self, text, color=None, attributes=None, link=None):
        """
        Initializes a table cell.
        :param text: Text
        :param color: Color
        :param attributes: Attributes
        :param link: Link to add to the cell text
        """
        self._as_text = text
        if color:
            color_attribute = [('class', color)]
            if attributes:
                attributes += color_attribute
            else:
                attributes = color_attribute
        if link:
            super(HtmlTableCell, self).__init__('td', None, attributes)
            self.add_html_object(HtmlElement('a', text, [('href', link)]))
        else:
            super(HtmlTableCell, self).__init__('td', text, attributes)

    @property
    def text(self):
        """
        Returns the text belonging to this tag.
        :return: Text
        """
        return self._as_text
