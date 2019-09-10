import logging
from typing import List, Tuple, Dict

import os

from camel.app.components.validation.typinghit import TypingHit


class MlstTabularParser(object):
    """
    This class is used to parse tabular output from the MLST tool.
    """

    @staticmethod
    def parse_tabular_input(file_path: str, include_imperfect: bool = False) -> List[Tuple[str, str]]:
        """
        Parses a tabular input file.
        :param file_path: File path
        :param include_imperfect: If True, non-perfect hits are included
        :return: List of (locus, allele_id)
        """
        with open(file_path) as handle:
            header = handle.readline()
            if header.split('\t')[0] != 'Locus' or header.split('\t')[1] != 'Allele':
                raise ValueError("Invalid tabular file: {}".format(file_path))
            allele_ids = [(hit.locus, hit.allele_id if hit.is_perfect() or include_imperfect else '-') for
                          hit in [TypingHit.parse(l, '\t') for l in handle.readlines()]]
            logging.debug('Nb. of alleles: {}'.format(len(allele_ids)))
        return allele_ids

    @staticmethod
    def parse_tabular_all(tabular_input_files: List[Tuple[str, str]], include_imperfect: bool)\
            -> Dict[str, List[Tuple[str, str]]]:
        """
        Parses a list of tabular input files.
        :param tabular_input_files: List of input files + file name
        :param include_imperfect: If True, imperfect hits are included
        :return: Dictionary of detected alleles by sample name
        """
        allele_ids = {}
        for tabular_file, file_name in tabular_input_files:
            logging.debug('Parsing file: {}'.format(file_name))
            try:
                sample_name = os.path.splitext(os.path.basename(file_name))[0]
                logging.debug('Sample: {}'.format(sample_name))
            except IndexError:
                raise ValueError("Cannot determine sample name from: {}".format(file_name))
            allele_ids[sample_name] = MlstTabularParser.parse_tabular_input(tabular_file, include_imperfect)
        return allele_ids
