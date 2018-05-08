import abc
import logging
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class DatabaseParser(object, metaclass=abc.ABCMeta):
    """
    Superclass for parsing FASTA databases.
    """

    def __init__(self):
        """
        Initializes a ResFinder parser.
        """
        self._records = []
        self.parsed_alleles = []

    def _create_generic_fasta(self, filename):
        """
        Creates a generic FASTA with seq_XX headers.
        :return: Absolute path to the generated FASTA file
        """
        new_records = []
        for i in range(0, len(self._records)):
            new_record = SeqRecord(self._records[i].seq, id='seq_{}'.format(str(i)), description='')
            new_records.append(new_record)
        with open(filename, 'w') as output_handle:
            SeqIO.write(new_records, output_handle, 'fasta')
        return os.path.abspath(filename)

    @abc.abstractmethod
    def parse(self, filename):
        """
        Parses the FASTA file.
        :param filename: FASTA file
        :return: Alleles
        """
        return

    @staticmethod
    @abc.abstractmethod
    def _get_annotation(allele_data):
        """
        Returns the annotation for the given allele.
        :param allele_data: Allele data
        :return: Annotation
        """
        return

    def _get_unique_full_id(self, allele_data):
        """
        Returns a identifier allele + seqid that is unique.
        :return: Full allele identifier
        """
        full_id = '{}_{}'.format(allele_data['allele_name'], allele_data['seq_id'])
        if full_id in self.parsed_alleles:
            logging.info('Duplicate allele: {}'.format(allele_data))
            allele_data['seq_id'] = self._get_unique_allele_seqid(allele_data['allele_name'])
            full_id = '{}_{}'.format(allele_data['allele_name'], allele_data['seq_id'])
        return full_id

    def _get_unique_allele_seqid(self, allele_name):
        """
        Adds a seq id to the allele name so the combination of allele name + seq_id is unique.
        :param allele_name: Allele name
        :return: None
        """
        seq_id = 1
        combination = '{}_{}'.format(allele_name, seq_id)
        new_seq_id = int(seq_id)
        while combination in self.parsed_alleles:
            new_seq_id += 1
            combination = '{}_{}'.format(allele_name, new_seq_id)
        return str(new_seq_id)
