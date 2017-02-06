from app.components.blasthit import DEFAULT_COLUMNS, SEQEXTRACTION_COLUMNS_WITH_SEQS, SEQ_COLUMNS
from app.components.blasthit.blastntsvhit import BlastnTSVHit
from app.components.blasthit.blastntsvhitwithseq import BlastnTSVHitWithSeq


class BlastnTSVFile(object):
    """
    Class to handle blastn tsv output file
    """
    # TODO remove funcs for backward compatibility

    def __init__(self, blastn_tsv, with_seq=False, columns=None):
        """
        Initialize the class
        """
        self.blastn_tsv = blastn_tsv
        self.with_seq = with_seq
        if self.with_seq:
            self.hit_class = BlastnTSVHitWithSeq
        else:
            self.hit_class = BlastnTSVHit
        if columns is None:
            if with_seq:
                self._columns = SEQEXTRACTION_COLUMNS_WITH_SEQS
            else:
                self._columns = DEFAULT_COLUMNS
        else:
            self._columns = columns

    def get_inform_column_names(self):
        # backward compatibility, will be removed
        return self.inform_columns()

    def get_column_names(self):
        # backward compatibility, will be removed
        return self.columns()

    @property
    def inform_columns(self):
        """
        Returns columns without SEQ_COLUMNS (e.g., for html output)
        :return: columns without SEQ_COLUMNS
        """
        if self.with_seq:
            columns = []
            columns += self.columns
            for col in SEQ_COLUMNS:
                columns.remove(col)
            return columns
        else:
            return self.columns

    @property
    def columns(self):
        """
        Returns column names
        :return: _columns
        """
        return self._columns

    def read_hits(self, key='qseqid'):
        """
        Returns hits information as hash
        :return: hits grouped by key
        """
        hits = {}
        with open(self.blastn_tsv, 'r') as hits_file:
            for hit_inform in hits_file.readlines():
                hit = self.hit_class(hit_inform, self._columns)
                seqid = hit.qseqid if key == 'qseqid' else hit.sseqid
                if seqid in hits:
                    hits[seqid].append(hit)
                else:
                    hits[seqid] = [hit]
        return hits

    def read_hits_as_list(self):
        """
        Returns hits information in a list
        :return: hits in a list
        """
        hits = []
        with open(self.blastn_tsv, 'r') as hits_file:
            for hit_inform in hits_file.readlines():
                hits.append(self.hit_class(hit_inform, self._columns))
        return hits
