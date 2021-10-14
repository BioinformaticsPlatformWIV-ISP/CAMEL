import logging
import random
from abc import ABC, abstractmethod
from typing import List, Dict, Any

from camel.app.components.blasthit.blastnhit import BlastnHit


class SegmentTyping(ABC):

    """
    Class for segment typing
    """

    def __init__(self, segment: str, blast_hits: List[BlastnHit], seed: int = None):
        """
        Initializes the object
        :param segment: Segment that the hits map to (can also be single_segment)
        :param blast_hits: List of Blast hits for the given segment
        :param seed: Optional seed used for tie breaking
        """
        self._segment = segment
        self._blast_hits = blast_hits
        self._best_candidate_targets = []
        self._seed = seed if seed else random.randint(1, 1000000)
        logging.info(f'No seed given for SegmentTyping, using {self._seed}')
        self._ambiguous = None

    @property
    @abstractmethod
    def stats(self) -> Dict[str, Any]:
        """
        Return statistics about the typing
        :return: Dictionary with statistics
        """
        pass

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
        """
        Returns the list of best candidate targets that were identified.
        :return: List of best candidate targets
        """
        return self._best_candidate_targets

    @property
    def ambiguous(self) -> bool:
        """
        Returns whether the results for the best hits were ambiguous.
        :return: Ambiguous results or not
        """
        return self._ambiguous

    @abstractmethod
    def _quality_check(self) -> None:
        """
        Quality check the typing results, log the possible issues, and set the ambiguous variable.
        - Low top hit count
        - Multiple best target sequences
        :return: None
        """
        pass
