from typing import List, Union

from camel.app.components.blasthit.blastnasnparser import BlastnAsnParser
from camel.app.components.blasthit.blastnhit import BlastnHit
from camel.app.components.seqid.seqidparser import SeqIDParser


class InfluenzaBlastnAsnParser(BlastnAsnParser):

    """
    Class that deals with parsing BLAST results from an ASN file specifically
    for Influenza.
    """

    def __init__(self, blastn_file: str, multi_segment: bool, seqid_parser_type: str, subtyping_method: str,
                 exclude_tax_columns: bool = True, folder: str = None):
        """
        Initializes the object
        :param blastn_file: File location of the BLASTn output file
        :param multi_segment: Is it a multi-segment genome or not
        :param seqid_parser_type: Type of the SeqID parser
        :param subtyping_method: Subtyping method used (i.e. blast or asm)
        :param exclude_tax_columns: Should columns that require a taxonomy db be excluded
        :param folder: Optional folder to run the commands in
        """
        super(InfluenzaBlastnAsnParser, self).__init__(blastn_file, exclude_tax_columns=exclude_tax_columns, folder=folder)
        self._check_subtyping_method(subtyping_method)
        self._multi_segment = multi_segment
        self._seqid_parser_type = seqid_parser_type
        self._refseq_id = 'sseqid' if subtyping_method == 'alignment' else 'qseqid'
        self._segment_hits = {}

    @staticmethod
    def _check_subtyping_method(subtyping_method):
        if subtyping_method not in {'assembly', 'alignment'}:
            raise ValueError(f'Given subtyping method is invalid for InfluenzaBlastnAsnParser: {subtyping_method}')

    def get_segment_hits(self, segment: str) -> Union[None, List[BlastnHit]]:
        """
        Returns all BLAST hits for the requested segment.
        :param segment: Segment to return BLAST hits for
        :return: List of BLASTn hits or None if no hits for the segment were found
        """
        return self._segment_hits[segment] if segment in self._segment_hits else None

    def group_hits_per_segment(self) -> None:
        """
        Groups all BLASTn hits per segment. If it is not a multi-segment genome,
        the hits are assigned to a 'dummy' segment called single_segment.
        :return: None
        """
        if self._multi_segment:
            for hit in self._hits:
                segment = SeqIDParser(getattr(hit, self._refseq_id), self._seqid_parser_type).segment
                if segment in self._segment_hits:
                    self._segment_hits[segment].append(hit)
                else:
                    self._segment_hits[segment] = [hit]
        else:
            self._segment_hits['single_segment'] = self._hits
