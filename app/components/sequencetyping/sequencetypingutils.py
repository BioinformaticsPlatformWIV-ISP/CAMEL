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
        :param fasta_file: FASTA file
        :return: None
        """
        with open(fasta_file) as handle:
            for seq in SeqIO.parse(handle, 'fasta'):
                m = re.match('.*([-_])\d+$', seq.id)
                if m is None:
                    raise ValueError("Cannot determine allele delimiter")
                delimiter = m.group(1)
                logging.info("Detected delimited: '{}'".format(delimiter))
                return delimiter
