from app.components.blasthit.blastntsvhit import BlastnTSVHit


class BlastnTSVHitWithSeq(BlastnTSVHit):

    """
    Class to handle customized blastn tabular output with algined sequences

    Supported blastn outfmt 6 columns:
    'qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send qseq sseq evalue bitscore strand qcovs qcovhsp'
    """
    # TODO MOVED TO PIPELINE SPECIFIC CODE
    SEQEXTRACTION_COLUMNS_WITH_SEQS = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "gaps", "qstart", "qend",
        "sstart", "send", "qseq", "sseq", "evalue", "bitscore", "strand", "qcovs", "qcovhsp"
    ]

    SEQ_COLUMNS = ['sseq', 'qseq']

    def __init__(self, hit, columns=BlastnTSVHitWithSeq.SEQEXTRACTION_COLUMNS_WITH_SEQS):
        """
        Initialize
        :param hit: the txt string describe hit
        """
        if any(col not in columns for col in BlastnTSVHitWithSeq.SEQ_COLUMNS):
            raise ValueError("BlastnTSVHitWithSeq requires both sseq and qseq to be reported.")

        super(BlastnTSVHitWithSeq, self).__init__(hit, columns)

        self._hit_str_without_seqs = "\t".join(self.__extract_none_seq_informs())

    def __extract_none_seq_informs(self):
        """
        Extracts hit inform without sequence fields (sseq, qseq)
        :return: seq informs in list
        """
        informs = []
        for x in self.columns:
            if x not in SEQ_COLUMNS:
                informs += [str(getattr(self, x))]
        return informs

    def print_with_seqs(self):
        """
        Print hit information with sequences
        """
        return self._hit_str

    def __str__(self):
        return self._hit_str_without_seqs

    def __repr__(self):
        return str(self)
