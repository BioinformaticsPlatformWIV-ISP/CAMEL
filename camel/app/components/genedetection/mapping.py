from typing import Dict

import json


class Mapping(object):
    """
    This class contains a mapping of the converted sequence name to the original header.
    The main purpose of this class is to avoid cluttering of the log with the complete mapping as a dictionary.
    """

    def __init__(self, content: Dict[str, str]):
        """
        Initializes a mapping.
        :param content: Mapping content
        """
        self._content = content

    @staticmethod
    def parse(input_file: str) -> 'Mapping':
        """
        Parses a mapping from a file.
        :param input_file: Input file
        :return: Mapping
        """
        with open(input_file) as handle:
            return Mapping(json.load(handle))

    @property
    def content(self) -> Dict[str, str]:
        """
        Returns the mapping content.
        :return: Content
        """
        return self._content

    def __repr__(self) -> str:
        """
        Returns the printable representation of the mapping.
        :return: Representation
        """
        return 'Mapping({} items)'.format(len(self._content))

    def get(self, key) -> str:
        """
        Returns the item with the given key.
        :param key: Key
        :return: Item
        """
        return self._content[key]
