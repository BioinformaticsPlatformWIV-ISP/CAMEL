import logging
import pprint
from typing import Dict, List, Union

from Bio.SeqIO import SeqRecord

from camel.app.camel import Camel
from camel.app.components.blasthit.blastnasnparser import BlastnAsnParser
from camel.app.components.blasthit.blastnhit import BlastnHit
from camel.app.components.blasthit.blastnhitindelscanner import BlastnHitIndelScanner
from camel.app.components.files.fastautils import FastaUtils
from camel.app.components.seqid.seqidparser import SeqIDParser
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.tools.tool import Tool


class BlastnSequenceInformExtractor(Tool):

    """
    Calculate consensus sequence statistics from Blastn ASN output

    INPUTs:
    - FASTA_REF: reference sequence(s) fasta file (for refseq length inform)
    - ASN: ASN Blastn output
    - SeqIDParser [optional]: for off-target hit checkup

    INFORMs (output):
    - Target_consensus_seq_inform: consensus sequence statistics

    Application Notes:
    - optional input SeqIDParser is set only for calculating stats for conseq.
    - optional input SeqIDParser does NOT applicable for calculating stats for inconsistent segment (contig id w/o segment information).
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialization
        :param camel: Camel instance
        :return: None
        """
        super().__init__('BlastnSequenceInformExtractor', '0.1', camel)
        self._check_offtarget_hit = False

    def _execute_tool(self) -> None:
        """
        Function to run the tool
        :return: None
        """
        self.__gather_sequence_statistics()

    def _check_input(self) -> None:
        """
        Check input specs
        :return: None
        """
        if 'ASN' not in self._tool_inputs:
            raise InvalidInputSpecificationError("BlastnSequenceInformExtractor required blastn ASN input file missing.")

        if 'FASTA_REF' not in self._tool_inputs:
            raise InvalidInputSpecificationError("BlastnSequenceInformExtractor required FASTA_REF input file (reference sequence fasta file) missing.")

        if 'seqIDParser_type' in self._parameters:
            self._check_offtarget_hit = True
            logging.info("Input SeqIDParser specified, off target hits detection on.")
        else:
            logging.info("No SeqIDParser input specified, off target hits detection off.")

        super()._check_input()

    def __gather_sequence_statistics(self) -> None:
        """
        Gather sequence statistics
        :return: None
        """
        self._informs['target_consensus_seq_inform'] = {}
        reference_seqs = FastaUtils.read_as_dict(self._tool_inputs['FASTA_REF'][0].path)
        blastn_hits = self.__retrieve_blastn_hits()

        for qseqid, hits in blastn_hits.items():
            logging.debug(f'Process query {qseqid} --------- ')

            # clean off-target hits if possible
            if self._check_offtarget_hit:
                hits = self.__filter_offtarget_hits(hits)
                if len(hits) == 0:
                    logging.debug(f'  No hits remain after filtering off-target hits for query {qseqid}.')
                    continue

            hit_informs = []
            query_rcd = self.__get_sequence_record(qseqid, reference_seqs)
            if len(hits) > 1:
                # NOTE: this can happen when there are no-coverage section in the sequences (which is represented as
                #       continuous 'N's). Then Blast will break it into different HSPs (local alignment!)
                logging.info(f'   WARNING: more than one hit for query {qseqid}.')
            logging.debug(f'  query hits: \n{pprint.pformat(hits, indent=2)}\n')
            for hit in hits:
                hit_informs.append(BlastnSequenceInformExtractor.__extract_hit_statistics(hit, len(query_rcd.seq)))
            if len(hits) > 0:
                seqid = hits[-1].sseqid if self._parameters['sequence_id'].value == 'sseqid' else hits[-1].qseqid
                self._informs['target_consensus_seq_inform'][seqid] = hit_informs

    def __retrieve_blastn_hits(self) -> Dict[str, List[BlastnHit]]:
        """
        Retrieve blastn hits and grouped by query (qseqid)
        :return: hits (as SequenceExtractionBlastnHitWithSeqs object) grouped by query
        """
        blastn_parser = BlastnAsnParser(self._tool_inputs['ASN'][0].path, seq_columns=True, exclude_tax_columns=True, folder=self._folder)
        hits_by_query = {}
        for hit in blastn_parser.hits:
            if hit.qseqid in hits_by_query:
                hits_by_query[hit.qseqid].append(hit)
            else:
                hits_by_query[hit.qseqid] = [hit]
        return hits_by_query

    @staticmethod
    def __get_sequence_record(seq_id: str, sequences: Dict[str, SeqRecord]) -> SeqRecord:
        """
        Retrieve sequence record of specific seq_id
        :param seq_id: Sequence id
        :param sequences: List of sequence records to search
        :return: SeqRecord of id
        """
        if seq_id not in sequences:
            raise KeyError(f'Could NOT identify the sequence with id {seq_id}.')

        return sequences[seq_id]

    def __filter_offtarget_hits(self, hits: List[BlastnHit]) -> List[BlastnHit]:
        """
        Check the segment information of the hit query and target and remove hits with different segments
        :param: Hits from same query segment
        :return: List of hits with off-target hits removed
        """
        filtered_hits = []
        cleaned_hits = []
        for hit in hits:
            if self.__is_offtarget_hit(hit):
                filtered_hits.append(hit)
            else:
                cleaned_hits.append(hit)

        if len(filtered_hits) > 0:
            logging.debug(f'  Off-target query hits (filtered): \n{pprint.pformat(filtered_hits, indent=2)}\n')
        else:
            logging.debug('  No off-target query hit found.')

        return cleaned_hits

    def __is_offtarget_hit(self, hit: BlastnHit) -> bool:
        """
        Check whether a hit is off-target (query and target are of different segments).

        Two types of sseqid possible:
        - assembly consensus sequence: A-Nicaragua-6585_04-2014-H3N2(CY240484)-MP-extracted
        - alignment consensus sequence: A-Passo_Fundo-LACENRS-1854-2015-H3N2(KY935475)-HA:1-1762

        :return: True if off-target
        """
        # Cleanup the region information if exist
        if ':' in hit.sseqid:
            sseqid = hit.sseqid.split(":")[0]
        elif '-extracted' in hit.sseqid:
            sseqid = hit.sseqid.replace('-extracted', '')
        else:
            sseqid = hit.sseqid

        query_segment = SeqIDParser(hit.qseqid, self._parameters['seqIDParser_type']).segment
        subject_segment = SeqIDParser(sseqid, self._parameters['seqIDParser_type']).segment

        return query_segment != subject_segment

    @staticmethod
    def __extract_hit_statistics(hit: BlastnHit, query_len: int) -> Dict[str, Union[str, int]]:
        """
        Extract hit statistics based on blastn hit.
        :param hit: blastn hit
        :param query_len: query sequence length
        :return: hit_inform dictionary contains hit statistics
        """
        hit_inform = {'coverage': hit.qcovs, 'start': hit.qstart, 'end': hit.qend, 'length': hit.length}
        # Note that 'coverage' statistics should remain the same across different hits if exists, as it is qcovs (considering
        # all hits) instead of hit specific qcovhsp.
        n_count = hit.sseq.count('N')
        # SNP count excludes 'N' bases
        hit_inform['snp_count'] = hit.mismatch - n_count
        hit_inform['n_count'] = n_count
        hit_inform['query_length'] = query_len
        insertions = BlastnHitIndelScanner.scan_sequence_indels(hit.qseq, '+', hit.qstart)
        if insertions:
            hit_inform['insertion_count'] = len(insertions)
            hit_inform['inserted_base_count'] = sum([x.length for x in insertions])
        else:
            hit_inform['insertion_count'] = 0
            hit_inform['inserted_base_count'] = 0
        deletions = BlastnHitIndelScanner.scan_sequence_indels(hit.sseq, '-', hit.qstart)
        if deletions:
            hit_inform['deletion_count'] = len(deletions)
            hit_inform['deleted_base_count'] = sum([abs(x.length) for x in deletions])
        else:
            hit_inform['deletion_count'] = 0
            hit_inform['deleted_base_count'] = 0
        return hit_inform
