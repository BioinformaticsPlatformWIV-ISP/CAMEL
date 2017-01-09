from app.components.html.element import HtmlElement


class HtmlTag(HtmlElement):
    """
    Class that represents a HTML element.
    """

    def __init__(self, tag_name, text, attributes):
        """
        Initialize a new HTML element.
        :param tag_name: HTML tag name
        :param text: Text displayed in the element
        :param attributes: Attributes
        """
        self._text_value = text
        self._tag_name = tag_name
        self._attributes = attributes
        super(HtmlTag, self).__init__()

    @property
    def text_value(self):
        """
        Text value of this tag.
        :return: Text value
        """
        return self._text_value

    @property
    def tag_name(self):
        """
        Returns the tag name
        :return: Tag name
        """
        return self._tag_name

    def _generate_html(self):
        """
        Generates the HTML representation of this element.
        :return: HTML representation
        """

        with self._get_tag(self._tag_name, self._attributes):
            self._text(self._text_value)
        return self._doc.getvalue()
