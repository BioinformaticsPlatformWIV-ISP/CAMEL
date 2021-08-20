import logging
from collections import Counter
from typing import List, Dict, Any

from camel.app.components.blasthit.blastnhit import BlastnHit
from camel.app.components.segmenttyping.segmenttyping import SegmentTyping


class SegmentTypingReads(SegmentTyping):

    """
    Class for segment typing based on blasting of reads against a reference database
    """

    def __init__(self, segment: str, blast_hits: List[BlastnHit], seed: int = None):
        """
        Initializes the object
        :param segment: Segment that the hits map to (can also be single_segment)
        :param blast_hits: List of Blast hits for the given segment
        :param seed: Optional seed used for tie breaking
        """
        super().__init__(segment, blast_hits, seed)
        self._max_cnt = None
        self._counts = None
        self._count_target_hits()
        self._quality_check()

    @property
    def stats(self) -> Dict[str, Any]:
        """
        Return statistics about the typing.
        :return: Dictionary with statistics
        """
        return {'refseqid': self.best_target,
                'candidates': self.best_candidate_targets,
                'ambiguous': self.ambiguous,
                'counts': self.counts}

    def counts(self):
        """
        Returns the count of the most commonly found sseqid.
        :return: Count of the most commonly found sseqid
        """
        return self._counts

    def _count_target_hits(self) -> None:
        """
        Counts the number of times each sseqid is identified in the blast hits. The most common
        ones are assigned to a variable as well as the maximum count that was observed.
        :return: None
        """
        cnt = Counter([h.sseqid for h in self._blast_hits])
        self._counts = cnt.most_common()
        self._max_cnt = max([c[1] for c in self._counts])
        self._best_candidate_targets = [c[0] for c in self._counts if c[1] == self._max_cnt]

    def _quality_check(self):
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