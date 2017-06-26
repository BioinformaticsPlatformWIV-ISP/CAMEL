import copy


from app.components.blasthit import BLASTN_INT_COLUMNS, BLASTN_FLOAT_COLUMNS, BLASTN_SEQ_COLUMNS
from app.error.invalidparametererror import InvalidParameterError


class BlastnFmt6TSVParser(object):

    """
    Class to parse blastn outfmt 6 generated TSV file and return hits (as hash or list)
    """
    DEFAULT_OUTFMT6_COLUMNS = ('qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                               'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')

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

    def __check_hit_valid_ity(self, hit_informs):
        """
        Check the validity of a blastn hit
        :param hit_informs: one blastn outfmt6 hit information in list
        :return: True if valid, raise ValueError otherwise
        """
        if len(hit_informs) != len(self._columns):
            raise ValueError(
                "The number of data of blastn hit differs from the number of headers specified. #Header {!r} vs #data {!r}. Headers {!r}, data {!r}".format(len(self._columns), len(hit_informs), self._columns, hit_informs))

    def __parse_blastn_hit_data(self, hit_data):
        """
        Parse blastn outfmt6 hit data (row) to create a hit (as dictionary)
        :param hit_data: blastn outfmt6 hit record (one row in tsv file)
        :return: blastn hit in dict
        """
        hit_informs = hit_data.split("\t")
        hit_informs[-1] = hit_informs[-1].rstrip()

        self.__check_hit_valid_ity(hit_informs)

        blastn_hit = {}
        for col_name, data in zip(self._columns, hit_informs):
            try:
                if col_name in BLASTN_INT_COLUMNS:
                    blastn_hit[col_name] = int(data)
                elif col_name in BLASTN_FLOAT_COLUMNS:
                    blastn_hit[col_name] = float(data)
                else:
                    blastn_hit[col_name] = data

            except ValueError:
                logging.error("INCOMPATIBLEDATA, Blastn hit data incompatible with its column specification, column name {!r} value {!r}.".format(
                    col_name, data))
                raise ValueError

        return blastn_hit

    def read_hits_as_hash(self, key='qseqid'):
        """
        Returns hits information as hash, grouped by 'key'.
        :param key: data column to group hits, qseqid/sseqid
        :return: blastn hits grouped by key in hash
        """
        if key not in ('qseqid', 'sseqid'):
            raise InvalidParameterError(
                "Function 'read_hits_as_hash' support only using qseqid/sseqid as key, {!r} is provided.".format(key))

        if key not in self._columns:
            raise InvalidParameterError(
                "Blastn output does not contain required key column {!r}. Columns available {!r}.".format(key, self._columns))

        hits = {}
        with open(self._blastn_tsv, 'r') as hits_file:
            for hit_data in hits_file:
                hit = self.__parse_blastn_hit_data(hit_data)
                seqid = hit[key]
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
            for hit_data in hits_file:
                hit = self.__parse_blastn_hit_data(hit_data)
        return hits
