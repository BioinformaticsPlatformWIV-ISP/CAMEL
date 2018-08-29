import logging
import re
from Bio import SeqIO


class SequenceTypingUtils(object):
    """
    This class contains utility functions for sequence typing.
    """

    @staticmethod
    def determine_delimiter(fasta_file):
        """
        Returns the delimiter that is used in the FASTA file. Supported delimiters are '-' and '_'.
        For example a FASTA file with the first entry '>abcZ_2' will return '_' as the delimiter.
        :param fasta_file: FASTA file
        :return: None
        """
        with open(fasta_file) as handle:
            for seq in SeqIO.parse(handle, 'fasta'):
                m = re.match('.*([-_])\d+$', seq.id)
                if m is None:
                    raise ValueError("Cannot determine allele delimiter")
                delimiter = m.group(1)
                logging.info("Detected delimiter: '{}'".format(delimiter))
                return delimiter

    @staticmethod
    def get_allele_id(complete_name: str, regex: str=None) -> str:
        """
        Returns the allele id from a complete name.
        :param complete_name: Complete name (e.g. abcZ_2)
        :param regex: Allele id regex (e.g. '\d+$')
        :return: Allele id
        """
        if regex is None:
            regex = '\d+$'
        m = re.findall(regex, complete_name)
        if not len(m) == 1:
            raise ValueError("Cannot determine allele identifier for '{}' (RE: {})".format(complete_name, regex))
        return m[0]
