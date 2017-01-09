import logging
import re

from Bio import SeqIO

from app.components.resistance.dbparser import DatabaseParser


class ResFinderParser(DatabaseParser):
    """
    Parses the ResFinder database.
    """
    def __init__(self):
        """
        Initializes a ResFinder parser.
        """
        super(ResFinderParser, self).__init__()

    def parse(self, filename):
        """
        Parses the FASTA file.
        :param filename: FASTA file
        :return: Alleles
        """
        alleles = []
        with open(filename, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                self._records.append(record)
                allele_data = ResFinderParser.parse_header(record.description)
                self.parsed_alleles.append(self._get_unique_full_id(allele_data))

                allele_data['record'] = record
                allele_data['annotation'] = ResFinderParser._get_annotation(allele_data)

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
        accession = ResFinderParser.__parse_annotation(header)
        header = header.replace('_{}'.format(accession), '')
        seq_id = ResFinderParser.__parse_seq_id(header)
        allele_name = header.replace('_{}'.format(seq_id), '')
        gene_name = ResFinderParser.__parse_gene_name(header)

        parsed_data = {
            'gene_name': gene_name,
            'allele_name': allele_name,
            'seq_id': seq_id,
            'accession': accession
        }
        return parsed_data

    @staticmethod
    def _get_annotation(allele_data):
        """
        Returns the annotation for the given allele.
        :param allele_data: Allele data
        :return: Annotation
        """
        return allele_data['accession']

    @staticmethod
    def __parse_annotation(header):
        """
        Extracts the annotation from a header.
        :param header: Header
        :return: Annotation
        """
        m = re.match('^.*?\d_((?:N[CZ]_)?[A-Z\d]+(?:\.\d)?)$', header)
        if not m:
            raise StandardError('Cannot retrieve annotation from: {}'.format(header))
        return m.group(1)

    @staticmethod
    def __parse_seq_id(header):
        """
        Extracts the seq id from a header.
        :param header: Header without annotation
        :return: Seq id
        """
        m = re.match('.*_(\d+)$', header)
        if not m:
            raise StandardError('Cannot retrieve seq id from: {}'.format(header))
        return m.group(1)

    @staticmethod
    def __parse_gene_name(header):
        """
        Extracts the gene name from a header.
        :param header: Header
        :return: Gene name
        """
        m = re.match('^([A-Za-z\d()\']+)(\d*)', header)
        if not m:
            raise StandardError('Cannot retrieve gene name from: {}'.format(header))
        return m.group(1)
