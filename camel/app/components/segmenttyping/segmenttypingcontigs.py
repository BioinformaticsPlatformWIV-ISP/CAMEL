import pprint
from typing import List, Dict, Any

from camel.app.components.blasthit.blastnhit import BlastnHit
from camel.app.components.segmenttyping.segmenttyping import SegmentTyping
from camel.app.loggers import logger


class SegmentTypingContigs(SegmentTyping):
    """
    Class for segment typing based on blasting of a reference database against contigs
    """

    def __init__(self, segment: str, blast_hits: List[BlastnHit], seed: int = None):
        """
        Initializes the object
        :param segment: Segment that the hits map to (can also be single_segment)
        :param blast_hits: List of Blast hits for the given segment
        :param seed: Optional seed used for tie breaking
        """
        super().__init__(segment, blast_hits, seed)
        self._target_hits = {}
        self._summarized_stats = {}
        self._normalized_stats = {}
        self._top20 = {}
        self._collect_target_hits()
        self._process_target_stats()
        self._select_top_20()
        self._detect_best_targets()

    @property
    def stats(self) -> Dict[str, Any]:
        """
        Return statistics about the typing.
        :return: Dictionary with statistics
        """
        return {'refseqid': self.best_target,
                'candidates': self.best_candidate_targets,
                'top20': self._top20,
                'ambiguous': self.ambiguous}

    def _quality_check(self) -> None:
        """
        Quality check the typing results and log the possible issues.
        - Multiple best target sequences
        :return: None
        """
        # multiple best targets
        if 1 < len(self.best_candidate_targets) < 4:
            logger.debug(f"Segment {self._segment} has more than one best target sequence ({self.best_candidate_targets}) with hits statistics {self._normalized_stats[self.best_target]}.")
            self._ambiguous = True
        elif len(self.best_candidate_targets) >= 4:
            logger.warning(f"Segment {self._segment} has too many best target sequences ({self.best_candidate_targets}) with hits statistics {self._normalized_stats[self.best_target]}.")
            self._ambiguous = True

    def _collect_target_hits(self):
        """
        Collect typing results grouped by target refseq (hit.qseqid) from blastn hits
        :return: None
        """
        for hit in self._blast_hits:
            if hit.qseqid in self._target_hits:
                self._target_hits[hit.qseqid].append(hit)
            else:
                self._target_hits[hit.qseqid] = [hit]
        logger.debug(f'Segment typing: blast hits collected for segment {self._segment}: \n{self._target_hits}')

    def _process_target_stats(self) -> None:
        """
        Processes and normalizes typing information per target refseq over all hits and writes
        the normalized values to the debug log.
        :return: None
        """
        for refseq in self._target_hits:
            self._summarize_hits(refseq)
            self._normalize_stats(refseq)

        self._output_hits_stats()

    def _summarize_hits(self, refseq: str) -> None:
        """
        Summarize information over hits.
        :param refseq: Refseq ID to summarize the stats for
        :return: None
        """
        self._summarized_stats[refseq] = {
            'total_mismatch': 0,
            'total_gapopen': 0,
            'total_gaps': 0,
            'total_qcovs': 0,
            'total_length': 0,
            'total_identical_bases': 0,
            'hsp_count': 0,
            'contig_count': 0,
            'qcoverage_coordinates': []}
        subjects_seen = set()
        self._summarized_stats[refseq]['lowest_q_start_pos'] = 999999
        self._summarized_stats[refseq]['highest_q_end_pos'] = 0
        self._summarized_stats[refseq]['qlen'] = self._target_hits[refseq][0].qlen
        for hit in self._target_hits[refseq]:
            self._summarized_stats[refseq]['total_mismatch'] += hit.mismatch
            self._summarized_stats[refseq]['total_gapopen'] += hit.gapopen
            self._summarized_stats[refseq]['total_gaps'] += hit.gaps
            if hit.sseqid not in subjects_seen:
                self._summarized_stats[refseq]['total_qcovs'] += hit.qcovs
                self._summarized_stats[refseq]['contig_count'] += 1
                subjects_seen.add(hit.sseqid)
            # NOTE do NOT calculate ID_bases from length and pident, accuracy issue due to two times of roundings
            # - 1st pident (is rounded)
            # - 2nd when calculate ID_bases to convert into int
            self._summarized_stats[refseq]['total_identical_bases'] += hit.length - hit.mismatch - hit.gaps
            self._summarized_stats[refseq]['hsp_count'] += 1
            self._calculate_qcoverage_coordinates(refseq, hit)
            if hit.qstart < self._summarized_stats[refseq]['lowest_q_start_pos']:
                self._summarized_stats[refseq]['lowest_q_start_pos'] = hit.qstart
            if hit.qend > self._summarized_stats[refseq]['highest_q_end_pos']:
                self._summarized_stats[refseq]['highest_q_end_pos'] = hit.qend
        self._calculate_total_alignment_length(refseq)
        # self._summarized_stats[refseq]['max_alignment_length'] = max(x.length for x in self._target_hits[refseq])
        # self._summarized_stats[refseq]['variant_rate'] = round((self._summarized_stats[refseq]['total_mismatch'] + 3 * self._summarized_stats[refseq]['total_gaps'])
        #                                                        / float(self._summarized_stats[refseq]['total_length']), 4)

    def _calculate_qcoverage_coordinates(self, refseq: str, hit: BlastnHit) -> None:
        """
        Calculates the query coordinates (= refseq coordinates) that are covered by the different subjects.
        :param refseq: Refseq id
        :param hit: Blast hit to check coverage coordinates for
        :return: None
        """
        overlaps = [(i[0], i[1]) for i in self._summarized_stats[refseq]['qcoverage_coordinates'] if hit.qstart < i[1] and hit.qend > i[0]]
        if len(overlaps) == 0:
            self._summarized_stats[refseq]['qcoverage_coordinates'].append((hit.qstart, hit.qend))
        else:
            start_coord = min(overlaps[0][0], hit.qstart)
            end_coord = max(overlaps[-1][1], hit.qend)
            for interval in overlaps:
                self._summarized_stats[refseq]['qcoverage_coordinates'].remove(interval)
            self._summarized_stats[refseq]['qcoverage_coordinates'].append((start_coord, end_coord))
        self._summarized_stats[refseq]['qcoverage_coordinates'].sort()

    def _calculate_total_alignment_length(self, refseq: str) -> None:
        """
        Calculates the total number of bases that are covered by the different subjects.
        :return: None
        """
        self._summarized_stats[refseq]['total_length'] = 0
        for interval in self._summarized_stats[refseq]['qcoverage_coordinates']:
            # +1 as the calculation uses positions
            self._summarized_stats[refseq]['total_length'] += (interval[1] - interval[0] + 1)

    def _normalize_stats(self, refseq: str) -> None:
        """
        Normalize statistics based on the total coverage (qcovs) for over-covered segments, as in this case the raw statistics can be misleading.
        Note that the normalized value is converted into float with 1 digit for counting statistics.
        This is a bit of a naive approach as in theory you could see that all contigs overlap at the same point
        with a very low number of variants at those positions while you are still dividing the total number of
        variants by the total coverage number. The inverse is also possible, i.e. all variants in the overlapping
        parts but you would only divide by a relatively small number as you look at the entire sequence. There is
        however no straightforward way to know where all the variants are located so this is a pragmatic solution.

        The normalization is only done for entries with qcov > 1 as qcov < 0 is dealt with when correcting the
        variant rate based on the alignment length.
        :param refseq: Refseq ID to normalize the stats for
        :return: None
        """
        self._normalized_stats[refseq] = {}
        for key, value in self._summarized_stats[refseq].items():
            if key in {'total_mismatch', 'total_gapopen', 'total_gaps', 'total_identical_bases'}:
                cov = self._summarized_stats[refseq]['total_qcovs']
                self._normalized_stats[refseq][key] = round(100.0 * value / cov, 1) if cov > 100 else value
            else:
                self._normalized_stats[refseq][key] = value

    def _output_hits_stats(self) -> None:
        """
        Nicely formatted normalized hit stats output
        :return: None
        """
        logger.debug("Segment refseq hits statistics (sum out of blastn hits):")
        stats_keys = sorted(self._normalized_stats[next(iter(self._normalized_stats))].keys())
        table_data = []
        for refseq, stats in self._normalized_stats.items():
            row_data = [refseq]
            for key in stats_keys:
                row_data.append(str(stats[key]))
            table_data.append(row_data)
        logger.debug(['Refseq ID'] + stats_keys)
        logger.debug(pprint.pformat(table_data, indent=2, width=160))

    def _calculate_corrected_counts(self, refseq: str, max_alignment_length: int) -> float:
        """
        Calculates a corrected mismatch and identical bases count. This means that it will add the number of bases from the query
        that were not aligned to the subject (i.e. not all bases aligned from the refseq while the contig is longer)
        as mismatches with a maximum of the query length. This last condition makes sure that a refseq (=query) is not
        penalized by just being shorter than another refseq.

        The variant rate is defined as (the total number of mismatches + 3 * the number of gap positions + 0.5 * the difference
        between the longest alignment and the query length) divided by the maximum alignment length.
        :param refseq: Refseq id to calculate the variant rate for
        :param max_alignment_length:
        :return:
        """
        qlen = self._normalized_stats[refseq]['qlen']
        max_aln = qlen if max_alignment_length > qlen else max_alignment_length
        # The total length is the total number of aligned bases => this calculates the gaps that exist
        additional_mismatches = max_aln - self._normalized_stats[refseq]['total_length']
        self._normalized_stats[refseq]['total_mismatch'] += additional_mismatches
        self._normalized_stats[refseq]['total_identical_bases'] -= additional_mismatches
        self._normalized_stats[refseq]['variant_rate'] = (self._normalized_stats[refseq]['total_mismatch'] +
                                                          3 * self._summarized_stats[refseq]['total_gaps'] +
                                                          0.5 * (max_alignment_length - qlen)) / float(max_alignment_length)

    def _detect_best_targets(self):
        """
        Goes through the 20 refseq hits with the lowest variant rate and selects the best one(s).
        :return: None
        """
        best_variant_rate = 1
        best_refseq = ''
        for refseq in self._top20.keys():
            variant_rate = self._normalized_stats[refseq]['variant_rate']
            if best_variant_rate > variant_rate:
                best_refseq = refseq
                best_variant_rate = variant_rate
                self._best_candidate_targets = [best_refseq]
            elif variant_rate == best_variant_rate:
                self._best_candidate_targets.append(refseq)

        logger.info(f"Best reference {self.best_candidate_targets}:")
        logger.info(f"Reference statistics: {self._normalized_stats[best_refseq]})")

    def _select_top_20(self) -> None:
        """
        Calculates the variant rates for the different hits and creates a top 20
        based on these values
        :return: None
        """
        variant_rate_refseq = []
        max_alignment_length = max(self._normalized_stats[refseq]['highest_q_end_pos'] - self._normalized_stats[refseq]['lowest_q_start_pos'] + 1
                                   for refseq in self._normalized_stats)
        for refseq, seq_stats in self._normalized_stats.items():
            print(f'{refseq}\t{max_alignment_length}')
            self._calculate_corrected_counts(refseq, max_alignment_length)
            variant_rate_refseq.append((seq_stats['variant_rate'], refseq))
        top20 = sorted(variant_rate_refseq)[:20]
        logger.debug(f"Top 20 targets: {top20}")
        for entry in top20:
            self._top20[entry[1]] = {'summarized_stats': self._summarized_stats[entry[1]],
                                     'normalized_stats': self._normalized_stats[entry[1]]}
