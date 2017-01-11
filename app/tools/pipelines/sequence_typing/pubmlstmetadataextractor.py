import logging

from app.components.pubmlst.pubmlstparser import PubMLSTParser
from app.tools.tool import Tool


class PubmlstMetadataExtractor(Tool):
    """
    This tool parses PubMLST allele pages for additional metadata.
    Reports the information from the 'linked data' on PubMLST in the informs with allele name as key.
    Example:
        {'penA_1': ['Penicillin range (penA_1)', '>0.06 - 1 (intermediate) [n=413];'],
         'rpoB_3': ['Rifampicin range (rpoB_9)', '<=1 (susceptible) [n=31];']}
    The genes that have to be processed can be given as parameter 'loci'.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super(PubmlstMetadataExtractor, self).__init__('Sequence Typing: PubMLST Metadata Extractor', '0.1', camel)

    def _execute_tool(self):
        """
        Runs this tool.
        :return: None
        """
        if 'hits' not in self._input_informs:
            logging.info("No best hits found in informs")
            return

        checked_loci = self._parameters['loci'].value.split(', ')
        for loci in self._input_informs['hits']:
            allele_info = self._input_informs['hits'][loci]
            if loci in checked_loci:
                self.__add_linked_data(loci, allele_info['allele_id'], allele_info['url'])

    def __add_linked_data(self, name, allele_id, allele_url):
        """
        Adds the linked data for this allele.
        :param name: Gene name
        :param allele_id: Allele identifier
        :param allele_url: Allele URL
        :return: None
        """
        allele_name = '_'.join([name, PubmlstMetadataExtractor.__clean_allele_id(allele_id)])
        try:
            label, value = PubMLSTParser.parse_page(allele_url)
            self.informs[allele_name] = ['{} ({})'.format(label, allele_name), value]
        except RuntimeError as err:
            logging.warning('Cannot retrieve linked data for {} ({}, {})'.format(name, err.message, allele_url))

    @staticmethod
    def __clean_allele_id(allele_id):
        """
        Removes mismatch (*) or uncertainty (?) indicators from the allele id.
        :param allele_id: Allele identifier
        :return: Allele id
        """
        return allele_id.replace('*', '').replace('?', '')
