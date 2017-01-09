import urllib
from bs4 import BeautifulSoup


class PubMLSTParser(object):
    """
    Parsers for PubMLST web pages.
    """

    @staticmethod
    def parse_page(url):
        """
        Retrieves the linked data from a PubMLST page.
        :param url: URL of the web page to parse
        :return: Linked data
        """
        try:
            response = urllib.urlopen(url).read()
        except:
            raise RuntimeError('URL does not exist')
        html = BeautifulSoup(PubMLSTParser.__cleanup_html(response))
        try:
            html_value = html.find(text='Linked data').findNext('dd')
            html_label = html.find(text='Linked data').findNext('dt')
            return PubMLSTParser.__cleanup_label(html_label.text), html_value.text.replace(' PubMLST isolates', '')
        except AttributeError:
            raise RuntimeError('Cannot find linked data')

    @staticmethod
    def __cleanup_html(html_code):
        """
        Cleans up HTML code.
        :param html_code: HTML code
        :return: Cleaned up HTML code
        """
        return html_code.replace('<dd><', '<dd>&lt;')

    @staticmethod
    def __cleanup_label(label):
        """
        Capitalizes the label and replaces the underscore by a space.
        :param label: Input label
        :return: Cleaned up label
        """
        return label.capitalize().replace('_', ' ')
