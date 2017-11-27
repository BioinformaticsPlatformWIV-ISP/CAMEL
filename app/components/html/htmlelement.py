from yattag import Doc

from app.components.html.htmlbase import HtmlBase


class HtmlElement(HtmlBase):
    """
    This class contains a custom HTML element.
    """

    def __init__(self, tag, text=None, attributes=None):
        """
        Initializes an HTML element.
        :param tag: Tag
        :param attributes: HTML attributes
        :param text: Text
        """
        super(HtmlElement, self).__init__()
        self._tag_name = tag
        self._attributes = attributes
        self._tag_text = text

    @property
    def text(self):
        """
        Returns the HTML element as plain text.
        :return: Text
        """
        return self._tag_text

    # noinspection PyArgumentList
    def to_html(self):
        """
        Converts this element to HTML code.
        :return: HTML code
        """
        if self._attributes is None:
            self._attributes = []
        doc, tag, text = Doc().tagtext()
        if len(self._doc.getvalue()) == 0 and self._tag_text is None:
            doc.stag(self._tag_name, *self._attributes)
        else:
            with tag(self._tag_name, *self._attributes):
                if self._tag_text is not None:
                    text(self._tag_text)
                if len(self._doc.getvalue()) != 0:
                    doc.asis(self._doc.getvalue())
        return doc.getvalue()
