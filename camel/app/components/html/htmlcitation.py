import json
from typing import Dict, Any

from camel.app.components.html.htmlbase import HtmlBase
from camel.resources import DIR_CITATIONS


class HtmlCitation(HtmlBase):
    """
    This class is used to format citations in HTML.
    :param citation_data: Citation data (parsed from RIS)
    :return: None
    """

    def __init__(self, citation_data: Dict[str, Any]) -> None:
        """
        Initializes the citation.
        :param citation_data: Citation data
        """
        super().__init__()
        with self.get_tag('div', attributes=[('class', 'citations')]):
            with self.get_tag('p'):
                # Volume + (number optional)
                journal_parts = [citation_data['alternate_title3'], f", {citation_data['volume']}"]
                if 'number' in citation_data:
                    journal_parts.append(f" ({citation_data['number']})")

                # Citation text
                self.add_text('. '.join([
                    '; '.join(citation_data['first_authors']),
                    citation_data['primary_title'],
                    f"In <i>{''.join(journal_parts)}</i>",
                    ''
                ]))
                # Citation DOI / link
                with self.get_tag('a', [('href', f"https://dx.doi.org/{citation_data['doi']}")]):
                    self.add_text(f"DOI: {citation_data['doi']}")

    @staticmethod
    def parse_from_json(json_basename: str) -> 'HtmlCitation':
        """
        Parses a HTML citation from a JSON file.
        :param json_basename: Basename for the JSON file
        :return: Citation
        """
        json_citation = DIR_CITATIONS / f'{json_basename}.json'
        with json_citation.open(encoding='utf-8') as handle:
            data = json.load(handle)
        return HtmlCitation(data)
