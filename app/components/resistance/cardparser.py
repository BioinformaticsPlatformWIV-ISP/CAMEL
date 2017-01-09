import logging
import re

from Bio import SeqIO

from app.components.resistance.dbparser import DatabaseParser


class CardParser(DatabaseParser):
    """
    Parses the CARD database.
    """

    def __init__(self):
        """
        Initializes a CARD parser.
        """
        super(CardParser, self).__init__()

    def parse(self, filename):
        """
        Parses the FASTA file.
        :param filename: FASTA file
        :return: None
        """
        alleles = []
        with open(filename, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                self._records.append(record)
                allele_data = CardParser.parse_header(record.description)

                self.parsed_alleles.append(self._get_unique_full_id(allele_data))

                allele_data['record'] = record
                allele_data['annotation'] = CardParser._get_annotation(allele_data)

                alleles.append(allele_data)
        logging.info('{} Alleles parsed'.format(len(alleles)))
        return alleles

    @staticmethod
    def parse_header(header):
        """
        Parses a header.
        :param header: Header
        :return: Allele
        """
        header = header.strip()
        if not len(header.split('|')) == 6:
            raise ValueError("Invalid number of '|' in: {}".format(header))

        header_parts = header.split('|')

        # Parse the gene name and species
        # m = re.match("^(.+) \[(.+)\]$", header_parts[5])
        # if not m:
        #     raise ValueError('Cannot parse: {}'.format(header.strip()))

        parsed_data = {
            'source': header_parts[0],
            'accession_gb': header_parts[1],
            'strand': header_parts[2],
            'location': header_parts[3],
            'accession_card': header_parts[4],
            'gene_name': header_parts[5],
            'allele_name': header_parts[5],
            'seq_id': '1'
        }
        return parsed_data

    @staticmethod
    def _get_annotation(allele_data):
        """
        Returns the annotation for the given allele.
        :param allele_data: Allele data
        :return: Annotation
        """
        return allele_data['accession_card']
