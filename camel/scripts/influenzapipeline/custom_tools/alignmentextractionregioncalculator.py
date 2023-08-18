import logging
from pathlib import Path
from typing import Dict, List, Tuple, TextIO

from camel.app.camel import Camel
from camel.app.components.blasthit.indel import Indel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class AlignmentSeqExtractionRegionCalculator(Tool):

    """
    Class that generates interval files containing reference genome regions to be extracted (for alignment-based
    sequence extraction only). It combines VCF indels information with segment gaps information (by samtools depth)
    and optional HaplotypeCaller BAM output information to define target regions.

    Input Informs:
    ---------------------
    - IndelScanner: from VCFIndelScanner, indels information
    - DepthStats: mapping statistics from ReadMapping (including gaps information)
    - HCBamout (optional): mapping statistics from HaplotypeCaller generated BAM (with local re-assembly)

    Outputs:
    ---------------------
    - TXT_intervals: GATK-style .internvals file containing all regions to be extracted using FastaAlternateReferenceMaker

    Output Informs:
    ---------------------
    - extracted_regions: regions per segment that will be extracted

    Note:
    - GATK interval file format: google 'GATK Collected FAQs about interval lists'
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize tool
        :param camel: Camel instance
        :return: none
        """
        super().__init__('AlignmentSeqExtractionRegionCalculator', '0.1', camel)
        self.indels = {}
        self.segment_gaps = {}
        self.nocov_regions = {}
        self._extraction_regions = {}
        self._refseq_len = {}

    def _execute_tool(self) -> None:
        """
        Function to run AlignmentSeqExtractionRegionCalculator
        :return: None
        """
        self._check_required_informs()
        self._retrieve_informs()
        # incorporate HC bamout information (local de novo assembly) if available
        if self._has_hc_bam:
            logger.info("HaplotypeCaller generated BAM file provided, re-assembled regions will be considered.")
            self.__add_hc_haplotype_regions()
        self._identify_nocoverage_regions()
        self._set_output()

    def _check_required_informs(self) -> None:
        """
        Check whether the required informs are present
        :return: None
        """
        if 'IndelScanner' not in self._input_informs:
            raise InvalidInputSpecificationError("AlignmentSeqExtractionRegionCalculator required 'IndelScanner' information is missing.")
        if 'DepthStats' not in self._input_informs:
            raise InvalidInputSpecificationError("AlignmentSeqExtractionRegionCalculator required 'DepthStats' information is missing.")

    def _retrieve_informs(self) -> None:
        """
        Retrieve required information from other steps
        :return: None
        """
        self.indels = self._input_informs['IndelScanner']['indels']
        self.segment_gaps = self._input_informs['DepthStats']['segment_gaps']
        self._refseq_len = self._input_informs['DepthStats']['refseq_length']

    @property
    def _has_hc_bam(self) -> bool:
        """
        Check whether information of the HaplotypeCaller generated BAM out (by 'HC -bamout' opt.) is available
        :return: boolean, True if the information is available.
        """
        return 'HCBamout' in self._input_informs and self._input_informs['HCBamout'] is not None

    def __add_hc_haplotype_regions(self) -> None:
        """
        Incorporate HC local de novo assembled Haplotype covered regions into segment_gaps by adapting gaps based on overlaps
        with those regions
        :return: None
        """
        logger.info("Adapt segment gaps based on HaplotypeCaller re-assembled regions (from HCBamout).")
        haplotype_cov_regions = self.__calculate_cov_regions_from_gaps(self._input_informs['HCBamout']['segment_gaps'])
        logger.info(f"de novo Haplotypes covered regions: {haplotype_cov_regions}")
        logger.info(f"origin segment_gaps: {self.segment_gaps}")

        for segment, gaps in self.segment_gaps.items():
            if len(gaps) == 0:
                # skip segment without gaps
                continue

            if segment in haplotype_cov_regions:
                logger.debug(f"Processing segment {segment} ------------- ")
                segment_gaps_padded = []
                for gap in gaps:
                    has_overlap = False
                    overlap_cov_regions = []

                    for cov in haplotype_cov_regions[segment]:
                        if self.__overlap_region(cov, gap):
                            overlap_cov_regions.append(cov)
                            has_overlap = True

                    if not has_overlap:
                        segment_gaps_padded.append(gap)
                    else:
                        logger.debug(f'  gap: {gap}, overlap Haplotype regions {overlap_cov_regions}')
                        new_gaps = self.__adapt_gap_with_cov_regions(overlap_cov_regions, gap)
                        logger.debug(f'  adapted gaps: {new_gaps}')
                        segment_gaps_padded += new_gaps

                self.segment_gaps[segment] = segment_gaps_padded

        logger.info(f"updated segment_gaps: {self.segment_gaps}")

    def __calculate_cov_regions_from_gaps(self, gaps: dict) -> Dict[str, List[Tuple[int, int]]]:
        """
        Calculate covered regions based on reported gaps (for HCBam).
        :param gaps: the gaps grouped per segment (dictionary)
        :return: covered regions grouped per segment
        """
        cov_regions = {}
        for segment, segment_gaps in gaps.items():
            cov_regions[segment] = []
            gap_end_pos = 0
            for gap in segment_gaps:
                if gap[0] == gap_end_pos + 1:
                    gap_end_pos = gap[1]
                elif gap[0] > gap_end_pos + 1:
                    cov_regions[segment].append((gap_end_pos + 1, gap[0] - 1))
                    gap_end_pos = gap[1]
                else:
                    logger.warning(
                        f"gap {gap} start before the current start pos {gap_end_pos} of cov region.")

            # last gap
            if gap_end_pos != 'end' and gap_end_pos + 1 < self._refseq_len[segment]:
                cov_regions[segment].append((gap_end_pos + 1, self._refseq_len[segment]))

        return cov_regions

    @staticmethod
    def __overlap_region(region1: Tuple[int, int], region2: Tuple[int, int]) -> bool:
        """
        Check whether two regions overlap with each other
        :param region1: First region to check
        :param region2: Second region to check
        :return: True if overlap, False otherwise
        """
        if region1[0] > region2[1]:
            return False
        elif region1[1] < region2[0]:
            return False

        return True

    @staticmethod
    def __adapt_gap_with_cov_regions(cov_regions: List[Tuple[int, int]], gap: Tuple[int, int]) -> List[Tuple[int, int]]:
        """
        Adapt the gap with all cov regions overlap with it.
        :param cov_regions: cov regions overlap with gap
        :param gap: the gap
        :return: new_gaps, the list of gaps adapted
        """
        # Note that a gap might overlap with multiple cov regions. Cov regions do not overlap each other, and are
        # ordered by start_pos. Hence, after each iteration, the last gap in new_gaps should be used in the next
        # iteration to check its overlap with next region.
        new_gaps = []
        for cov_region in cov_regions:
            if new_gaps:
                gap = new_gaps.pop()
            new_gaps += AlignmentSeqExtractionRegionCalculator.__adapt_gap(cov_region, gap)

        return new_gaps

    @staticmethod
    def __adapt_gap(cov: Tuple[int, int], gap: Tuple[int, int]) -> List[Tuple[int, int]]:
        """
        Adapt a gap according to the cov regions
        :param cov: the covered region that overlaps with gap
        :param gap: the gap need to be handled
        :return: a list of new gaps considering cov region
        """
        if gap[0] >= cov[0] and gap[1] <= cov[1]:
            # cov: ---------------
            # gap:     ------
            return []
        elif gap[0] < cov[0] and gap[1] <= cov[1]:
            # cov:     ----------------
            # gap:  ---------
            return [(gap[0], cov[0] - 1)]
        elif gap[0] >= cov[0] and gap[1] > cov[1]:
            # cov:  ---------
            # gap:        ---------
            return [(cov[1] + 1, gap[1])]
        elif gap[0] < cov[0] and gap[1] > cov[1]:
            # cov:     --------
            # gap:  ----------------
            return [(gap[0], cov[0] - 1), (cov[1] + 1, gap[1])]

    def _identify_nocoverage_regions(self) -> None:
        """
        Calculate the real no-coverage region considering known indels (VCF), and gaps (readmapping)
        :return: None
        """
        # segment_gaps: the gaps (no coverage bases) detected by samtools depth statistics per segment
        #       indels: the indels detected by variance caller
        #
        # NOTE that segment_gaps is 0 coverage places, it is often much more strict than indels
        #
        # Definition & Coordinates:
        #   gap: (start, end) closed definition, start and end inclusive. 1-based position coordinate (1st base start at 1)
        # indel: (type, pos_of_proceeding_base, length) 1-based
        logger.info("  -- \t start searching overlapping gaps and indels ---------------")
        for segment, gaps in self.segment_gaps.items():
            if segment in self.indels:
                for gap in gaps:
                    gap_start = int(gap[0])
                    gap_length = int(gap[1]) - gap_start + 1
                    for ind in self.indels[segment]:
                        if ind.overlap(Indel('-', gap_start - 1, gap_length)):
                            if ind.type == '+':
                                logger.info("WARNNING segment gap overlap with insertion")
                            logger.info(f" segment {segment} overlap: gap {gap}, indel {ind}")
        logger.info("  -- \t finish searching overlapping gaps and indels ---------------")

        # Generally three categories of overlapping between gap and indel:
        #
        # Category I: indels inside of gaps (wrong indels)
        # --------------------------------------------------
        # - overlap from GATK VCF for sample 12-3663 (salmonella Assia)
        # 2016-11-21 22:55:32,333 - app.loggers - DEBUG - overlap: gap (2854496, 2854562), indel (-, 2854500, -1)
        # 2016-11-21 22:55:32,334 - app.loggers - DEBUG - overlap: gap (2854496, 2854562), indel (+, 2854504, 7)
        # 2016-11-21 22:55:32,334 - app.loggers - DEBUG - overlap: gap (2854496, 2854562), indel (+, 2854507, 1)
        # 2016-11-21 22:55:32,334 - app.loggers - DEBUG - overlap: gap (2854496, 2854562), indel (+, 2854508, 1)
        # 2016-11-21 22:55:32,334 - app.loggers - DEBUG - overlap: gap (2854496, 2854562), indel (+, 2854523, 4)
        # 2016-11-21 22:55:32,334 - app.loggers - DEBUG - overlap: gap (2854496, 2854562), indel (+, 2854539, 2)
        # 2016-11-21 22:55:32,334 - app.loggers - DEBUG - overlap: gap (2854496, 2854562), indel (+, 2854546, 2)
        # 2016-11-21 22:55:32,334 - app.loggers - DEBUG - overlap: gap (2854496, 2854562), indel (-, 2854549, -2)
        # 2016-11-21 22:55:32,334 - app.loggers - DEBUG - overlap: gap (2854601, 2871174), indel (-, 2854600, -2)
        # 2016-11-21 22:55:32,334 - app.loggers - DEBUG - overlap: gap (2854601, 2871174), indel (+, 2854605, 2)
        #
        # Category II: gaps over one side of indels (NONE observed)
        # --------------------------------------------------
        #
        # Category III: indels over the gaps
        # --------------------------------------------------
        # - GATK: 2016-11-21 22:55:32,335 - app.loggers - DEBUG - overlap: gap (2885617, 2885669), indel (-, 2885597, -82)
        # - no overlap is observed between gaps and indels reported by samtools and/or freebayes

        # Based on the observations, following strategy is taken:
        #  - for Category I & II indels, gap is maintained and indels will be ignored as gap region will not be extracted
        #  - for Category III indels, only one observed, and samtools and freebayes does not call indel at the same
        #    location. For better consistency over different callers, the gap is maintained and hence indel is ignored
        #
        # In summary, all gaps will be kept and overlapping indels ignored
        self.nocov_regions = self.segment_gaps

    def __output_regions_file(self, outfile_name: Path) -> None:
        """
        Output file specifying regions to be extracted based on nocov_regions
        :param outfile_name: filename for the output file
        :return: None
        """
        # region: (interval) specification follows samtools-style intervals with 1-based position coordinates
        #
        # gap: (start, end) closed definition, start and end inclusive. 1-based
        #      position coordinate (1st base start at 1)
        with open(outfile_name, 'w') as outf:
            for seqid in self.nocov_regions.keys():
                self._extraction_regions[seqid] = []

                if len(self.nocov_regions[seqid]) == 0:
                    # If no nocov regions, extract the whole sequence
                    refseq_end = self._refseq_len[seqid]
                    self.__output_extraction_region(seqid, 1, refseq_end, outf)
                else:
                    # If has nocov regions, extract regions between gaps
                    self.__output_extraction_region_with_gaps(seqid, self.nocov_regions[seqid], outf)

    def __output_extraction_region(self, seqid: str, start: int, end: int, outf: TextIO) -> None:
        """
        Output sequence and the region to be extracted into output file
        :param seqid: sequence id
        :param start: the start of the region to be extracted
        :param end: the end of the region to be extracted
        :param outf: the output file
        :return: None
        """
        outf.write(f"{seqid}:{start}-{end}\n")
        self._extraction_regions[seqid].append((start, end))

    def __output_extraction_region_with_gaps(self, seqid: str, gaps: List[Tuple[int, int]], outf: TextIO):
        """
        Output the regions to be extracted. Note that an extracted region is located either before the first gap, after the
        last gap, or in between two gaps
        :param seqid: sequence id
        :param gaps: gaps in the sequence
        :param outf: the output file
        :return: None
        """
        gap_end = 0
        for gap_region in gaps:
            region_start = gap_end + 1
            gap_start = int(gap_region[0]) - 1
            if gap_start > region_start:
                self.__output_extraction_region(seqid, region_start, gap_start, outf)
            if gap_region[1] == 'end':
                logger.warning(f"open-end gap {gap_region!r} found, the tail sequence in gap is skipped.")
                gap_end = gap_region[1]
            else:
                gap_end = int(gap_region[1])

        # Handling the possible tail segment after the last gap:
        #
        # 1. check whether there is any bases between 'gap_end' and 'refseq_end'
        # 2. if so, add the last segment (from the end of the last gap till the end of the sequence)
        if gap_end != 'end':
            refseq_end = self._refseq_len[seqid]
            if gap_end + 1 <= refseq_end:
                self.__output_extraction_region(seqid, gap_end + 1, refseq_end, outf)
            elif gap_end > refseq_end:
                logger.warning(
                    "gap ends further than the end of the sequences: gap_end {}, refseq_end {}".format(gap_end, refseq_end))

    def _set_output(self) -> None:
        """
        Set the output specification
        :return: None
        """
        outfile_name = self.folder / self._parameters['output'].value
        self._tool_outputs['TXT_intervals'] = [ToolIOFile(outfile_name)]
        self.__output_regions_file(outfile_name)
        self._informs['extracted_regions'] = self._extraction_regions
