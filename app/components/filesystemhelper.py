import unicodedata

import re


class FileSystemHelper(object):
    @staticmethod
    def make_valid(value):
        """
        Converts arbitrary strings to URL- and filename friendly values.
        :param value: Input value
        :return: URL- and filename friendly value
        """
        value = unicodedata.normalize('NFKD', unicode(value)).encode('ascii', 'ignore')
        value = unicode(re.sub('[^\w\s-]', '', value).strip())
        value = unicode(re.sub('[-\s]+', '-', value))
        return value
