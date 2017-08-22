from app.components.blasthit import BLASTN_SEQ_COLUMNS
from app.components.blasthit.blastntsvhit import BlastnTSVHit


class BlastnTSVHitWithSeq(BlastnTSVHit):

    """
    Class to handle customized blastn hit containing algined hit sequence(s) (sseq/qseq), providing convenient output
    formating by skipping sequence(s) by default.
    """

    def __init__(self, hit, columns):
        """
        Initialize BlastnTSVHitWithSeq
        :param hit: the txt string describe hit
        """
        super(BlastnTSVHitWithSeq, self).__init__(hit, columns)

        self._hit_str_without_seqs = "\t".join(self.__extract_none_seq_informs())

    def _check_validity(self, hit_inform):
        """
        Check the validity of the blast hit & its columns. If fails, raise error.
        :param hit_inform: hit information in list
        """
        # the length of data column and headers should match
        if len(hit_inform) != len(self._columns):
            raise ValueError(
                "The number of data of blastn hit differs from the number of headers specified. Headers {!r}, data {!r}".format(self._columns, hit_inform))

        # should contain some sequence data column
        if not any(x in self._columns for x in BLASTN_SEQ_COLUMNS):
            raise ValueError(
                "BlastnTSVHitWithSeq class is specialized to handle blastn output with sequence information (sseq/qseq specifier). Otherwise, use BlastnTSVHit class instead.")

    def __extract_none_seq_informs(self):
        """
        Extracts hit inform without sequence fields (sseq, qseq).
        :return: seq informs in list
        """
        informs = []
        for x in self._columns:
            if x not in BLASTN_SEQ_COLUMNS:
                informs += [str(getattr(self, x))]
        return informs

    def print_with_seqs(self):
        """
        Print hit information with sequences.
        :return: complete hit information in a string with hit sequences
        """
        return self._hit_str

    def __str__(self):
        return self._hit_str_without_seqs
