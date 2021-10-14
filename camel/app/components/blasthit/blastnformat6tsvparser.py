import copy
from typing import List, Union, Dict

from camel.app.error.invalidparametererror import InvalidParameterError


class BlastnFormat6TSVParser(object):

    """
    Class to parse blastn output format 6 generated TSV file and return hits (as hash or list)
    """

    DEFAULT_OUTFMT6_COLUMNS = ('qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                               'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
    BLASTN_INT_COLUMNS = ('qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'score', 'length', 'nident', 'mismatch',
                          'positive', 'gapopen', 'gaps', 'qcovs', 'qcovhsp')
    BLASTN_FLOAT_COLUMNS = ('pident', 'ppos', 'bitscore')
    BLASTN_SEQ_COLUMNS = ('sseq', 'qseq')

    def __init__(self, blastn_tsv: str, columns: List[str] = None) -> None:
        """
        Initialize the BlastnFmt6TSVParser
        :param blastn_tsv: blastn tsv output file as input
        :param columns: data columns specified for outfmt 6. Default one is provided, for customized run, this should be provided.
        :return: None
        """
        self._blastn_tsv = blastn_tsv

        if columns is None:
            self._columns = BlastnFormat6TSVParser.DEFAULT_OUTFMT6_COLUMNS
        else:
            self._columns = columns

        self.with_seq = self.__has_seq_in_columns()

    @property
    def inform_columns(self) -> List[str]:
        """
        Returns columns without BLASTN_SEQ_COLUMNS (e.g., for html output) for convenient data output. To output
        sequences is not desirable in most cases as it is hard to have a nice format.
        :return: columns without BLASTN_SEQ_COLUMNS
        """
        if self.with_seq:
            columns = copy.copy(self._columns)
            for col in BlastnFormat6TSVParser.BLASTN_SEQ_COLUMNS:
                columns.remove(col)
            return columns
        else:
            return self._columns

    @property
    def columns(self) -> List[str]:
        """
        Returns column names
        :return: column names
        """
        return self._columns

    def __has_seq_in_columns(self) -> bool:
        """
        Check whether the blastn output contains hit sequence(s) information
        :return: boolean, True if contains hit sequence, False otherwise
        """
        return any(x in self._columns for x in BlastnFormat6TSVParser.BLASTN_SEQ_COLUMNS)

    def __check_hit_validity(self, hit_informs: List[Union[str, int, float]]) -> None:
        """
        Check the validity of a blastn hit by verifying that the number of requested columns
        is equal to the columns present in the data
        :param hit_informs: one blastn outfmt6 hit information in list
        :return: None if valid, raise ValueError otherwise
        """
        if len(hit_informs) != len(self._columns):
            raise ValueError(f"The number of data of blastn hit differs from the number of headers specified. "
                             f"#Header {len(self._columns)} vs #data {len(hit_informs)}. Headers {self._columns}, "
                             f"data {hit_informs}")

    def __parse_blastn_hit_data(self, hit_data: str) -> Dict[str, Union[str, int, float]]:
        """
        Parse blastn outfmt6 hit data (row) to create a hit (as dictionary)
        :param hit_data: blastn outfmt6 hit record (one row in tsv file)
        :return: blastn hit in dict
        """
        hit_informs = hit_data.split("\t")
        hit_informs[-1] = hit_informs[-1].rstrip()

        self.__check_hit_validity(hit_informs)

        blastn_hit = {}
        for col_name, data in zip(self._columns, hit_informs):
            try:
                if col_name in BlastnFormat6TSVParser.BLASTN_INT_COLUMNS:
                    blastn_hit[col_name] = int(data)
                elif col_name in BlastnFormat6TSVParser.BLASTN_FLOAT_COLUMNS:
                    blastn_hit[col_name] = float(data)
                else:
                    blastn_hit[col_name] = data
            except ValueError:
                raise ValueError(f"Incompatible data, Blastn hit data incompatible with its column specification, "
                                 f"column name {col_name} value {data}.")

        return blastn_hit

    def read_hits_as_hash(self, key: str = 'qseqid') -> Dict[Union[str, int, float], List[Dict[str, Union[str, int, float]]]]:
        """
        Returns hits information as hash, grouped by 'key'.
        :param key: data column to group hits, qseqid/sseqid
        :return: blastn hits grouped by key in hash
        """
        if key not in {'qseqid', 'sseqid'}:
            raise InvalidParameterError(
                f"Function 'read_hits_as_hash' supports only using qseqid/sseqid as key, {key} is provided.")

        if key not in self._columns:
            raise InvalidParameterError(
                f"Blastn output does not contain required key column {key}. Columns available {self._columns}.")

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

    def read_hits_as_list(self) -> List[Dict[str, Union[str, int, float]]]:
        """
        Returns hits information in a list.
        :return: list of blastn hits
        """
        hits = []
        with open(self._blastn_tsv, 'r') as hits_file:
            for hit_data in hits_file:
                hit = self.__parse_blastn_hit_data(hit_data)
                hits.append(hit)
        return hits
