import abc
from abc import ABCMeta

from yattag import Doc


class HtmlElement(object):
    """
    Base class for HTML tags.
    """

    __metaclass__ = ABCMeta

    def __init__(self):
        """
        Initializes this HTML tag.
        """
        self._doc, self._tag, self._text = Doc().tagtext()
        self._html = self._generate_html()

    @property
    @abc.abstractmethod
    def text_value(self):
        """
        Textual value of the HtmlElement.
        :return: Text value
        """
        pass

    @abc.abstractmethod
    def _generate_html(self):
        """
        Generates the HTML code for this tag.
        :return: HTML code
        """
        pass

    @property
    def html_code(self):
        """
        Returns the HTML code.
        :return: HTML code
        """
        return self._html

    def _get_tag(self, tag, attributes=None):
        """
        Returns a tag with the given attributes.
        :param tag: Tag name
        :param attributes: HTML attributes
        :return: Tag
        """
        if attributes:
            return self._tag(tag, *attributes)
        else:
            return self._tag(tag)
