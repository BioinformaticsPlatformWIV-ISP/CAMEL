import logging
from typing import List

from camel.app.components.blasthit.blastnhit import BlastnHit
from camel.app.components.blasthit.indel import Indel


class BlastnHitIndelScanner(object):

    """
    Class to scan a blastn hit (with seq) to extract indels
    """

    @staticmethod
    def scan_sequence_indels(sequence: str, indel_type: str, start: int) -> List[Indel]:
        """
        Scan a sequence to extract its indels
        :param sequence: the raw sequence output of blastn hit
        :param indel_type: character, '-' deletion, '+' insertion
        :param start: the start position of the sequence on query (blastn hit qstart)
        """
        indels = []
        current_pos = start
        start_pos = 0
        in_indel = False
        for c in sequence:
            if c == '-':
                if not in_indel:
                    in_indel = True
                    start_pos = current_pos
            elif in_indel:
                # NOTE in this case: c != '-'
                indels.append(Indel(indel_type, start_pos, current_pos - start_pos))
                in_indel = False
            current_pos += 1
        return indels

    @staticmethod
    def scan_indels(blastn_hit: BlastnHit) -> List[Indel]:
        """
        Gather all indels in a blastn hit
        :param blastn_hit: a BlastnTSVHitWithSeq object containing a blastn hit information
        :return: indels, list of indels of class Indel
        """
        indels = []
        if blastn_hit.gaps == 0:
            return indels
        else:
            if blastn_hit.sseq and blastn_hit.qseq:
                # insertion on query sequence
                indels = BlastnHitIndelScanner.scan_sequence_indels(blastn_hit.qseq, '+', blastn_hit.qstart)
                # deletion on subject sequence
                indels += BlastnHitIndelScanner.scan_sequence_indels(blastn_hit.sseq, '-', blastn_hit.qstart)
                return indels
            else:
                logging.error("Blastnhit without sequences cannot be scanned format indels.")
                raise ValueError('Blastnhit contain gaps without sequences. IndelScanner fails to extract indels.')
