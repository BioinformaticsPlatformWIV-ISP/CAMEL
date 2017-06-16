import copy


from app.components.blasthit import BLASTN_SEQ_COLUMNS
from app.components.blasthit.blastntsvhit import BlastnTSVHit
from app.components.blasthit.blastntsvhitwithseq import BlastnTSVHitWithSeq
from app.error.invalidparametererror import InvalidParameterError


class BlastnFmt6TSVParser(object):

    """
    Class to parse blastn outfmt 6 generated TSV file and return hits (as hash or list)
    """
    # TODO remove old funcs in favor of properties
    DEFAULT_OUTFMT6_COLUMNS = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split(
        " ")

    def __init__(self, blastn_tsv, columns=None):
        """
        Initialize the BlastnFmt6TSVParser
        :param blastn_tsv: blastn tsv output file as input
        :param columns: the data columns specified for blastn outfmt 6. Default one is provided, for customized run, this should be provided.
        :return: None
        """
        self._blastn_tsv = blastn_tsv

        if columns is None:
            self._columns = BlastnFmt6TSVParser.DEFAULT_OUTFMT6_COLUMNS
        else:
            self._columns = columns

        self.with_seq = self.__has_seq_in_columns()
        if self.with_seq:
            self.hit_class = BlastnTSVHitWithSeq
        else:
            self.hit_class = BlastnTSVHit

    @property
    def inform_columns(self):
        """
        Returns columns without BLASTN_SEQ_COLUMNS (e.g., for html output) for convenient data output. To output sequences is not
        desirable in most cases as it is hard to have a nice format.
        :return: columns without BLASTN_SEQ_COLUMNS
        """
        if self.with_seq:
            columns = copy.copy(self._columns)
            for col in BLASTN_SEQ_COLUMNS:
                columns.remove(col)
            return columns
        else:
            return self._columns

    @property
    def columns(self):
        """
        Returns column names
        :return: column names
        """
        return self._columns

    def __has_seq_in_columns(self):
        """
        Check whether the blastn output contains hit sequence(s) information
        :return: boolean, True if contains hit sequence, False otherwise
        """
        return any(x in self._columns for x in BLASTN_SEQ_COLUMNS)

    def read_hits_as_hash(self, key='qseqid'):
        """
        Returns hits information as hash, grouped by 'key'.
        :param key: data column to group hits, qseqid/sseqid
        :return: blastn hits grouped by key in hash
        """
        if key not in ('qseqid', 'sseqid'):
            raise Inva2lidParameterError(
                "Function 'read_hits_as_hash' support only using qseqid/sseqid as key, {} is provided.".format(key))

        hits = {}
        with open(self._blastn_tsv, 'r') as hits_file:
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
        Returns hits information in a list.
        :return: list of blastn hits
        """
        hits = []
        with open(self._blastn_tsv, 'r') as hits_file:
            for hit_inform in hits_file.readlines():
                hits.append(self.hit_class(hit_inform, self._columns))
        return hits

    # TODO remove old functions in favor of properties

    # def get_inform_column_names(self):
    # backward compatibility, will be removed
    #     return self.inform_columns()

    # def get_column_names(self):
    # backward compatibility, will be removed
    #     return self.columns()
