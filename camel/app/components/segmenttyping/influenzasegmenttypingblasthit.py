import logging
import random
from collections import Counter
from typing import List

from camel.app.components.blasthit.blastnhit import BlastnHit


class InfluenzaSegmentTypingBlasthit(object):

    """
    Class to analyze BLASTn hits for Influenza segment typing
    """

    def __init__(self, segment: str, blast_hits: List[BlastnHit], seed: int = None):
        self._segment = segment
        self._blast_hits = blast_hits
        self._best_candidate_targets = None
        self._seed = seed if seed else random.randint(1, 1000000)
        logging.info(f'No seed given for InfluenzaSubtypingBlastHit, using {self._seed}')
        self._ambiguous = None
        self._max_cnt = None
        self._counts = None
        self.__count_target_hits()
        self.__quality_check()

    @property
    def best_target(self) -> str:
        """
        Returns the best blast hit that was found. In case of a tie, a random one is returned
        but the seed is set to return the same one every time.
        :return: Sseqid of the best hit
        """
        if len(self._best_candidate_targets) == 1:
            return self._best_candidate_targets[0]
        else:
            random.seed(self._seed)
            return self._best_candidate_targets[random.randint(0, len(self._best_candidate_targets)-1)]

    @property
    def best_candidate_targets(self) -> List[str]:
        return self._best_candidate_targets

    @property
    def ambiguous(self):
        return self._ambiguous

    @property
    def counts(self):
        return self._counts

    def __count_target_hits(self) -> None:
        """
        Counts the number of times each sseqid is identified in the blast hits. The most common
        ones are assigned to a variable as well as the maximum count that was observed.
        :return: None
        """
        cnt = Counter([h.sseqid for h in self._blast_hits])
        self._counts = cnt.most_common()
        self._max_cnt = max([c[1] for c in self._counts])
        self._best_candidate_targets = [c[0] for c in self._counts if c[1] == self._max_cnt]

    def __quality_check(self):
        """
        Quality check the typing results, log the possible issues, and set the ambiguous variable.
        - Low top hit count
        - Multiple best target sequences
        :return: None
        """
        # low max hit
        if self._max_cnt < 10:
            logging.warning(f"Subtyping: number of top hits is less than 10, only {self._max_cnt} hits.")

        # multiple best targets
        if 1 < len(self._best_candidate_targets) < 6:
            logging.debug(f"Subtyping: segment {self._segment} has more than one target sequences "
                          f"({self._best_candidate_targets}) with the maximal hit count {self._max_cnt}.")
            self._ambiguous = True
        elif len(self._best_candidate_targets) >= 6:
            logging.warning(f"Subtyping: segment {self._segment} has too many target sequences ({self._best_candidate_targets}) with the maximal hit count {self._max_cnt}.")
            self._ambiguous = True
        else:
            self._ambiguous = False
