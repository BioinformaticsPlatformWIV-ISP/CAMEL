import logging

from app.components.blasthit import BLASTN_INT_COLUMNS, BLASTN_FLOAT_COLUMNS, BLASTN_SEQ_COLUMNS


class BlastnTSVHit(object):

    """
    Class to handle blastn outfmt 6 hit w/o sequence information (sseq/qseq).
    """

    def __init__(self, hit, columns):
        """
        Initialize BlastnTSVHit instance
        :param hit: the txt string describe hit
        :param columns: specified data column names
        :return: None
        """
        hit_inform = hit.split("\t")
        hit_inform[-1] = hit_inform[-1].rstrip()
        self._columns = columns
        self._hit_str = hit
        self._check_validity(hit_inform)
        for col_name, data in zip(columns, hit_inform):
            try:
                if col_name in BLASTN_INT_COLUMNS:
                    setattr(self, col_name, int(data))
                elif col_name in BLASTN_FLOAT_COLUMNS:
                    setattr(self, col_name, float(data))
                else:
                    setattr(self, col_name, data)

            except ValueError:
                logging.error("Blastn hit data incompatible with its column specification, column name {!r} value {!r}.".format(
                    col_name, data))
                raise

    def _check_validity(self, hit_inform):
        """
        Check the validity of the blast hit & its columns. If fails, raise error.
        :param hit_inform: hit information in list
        """
        # the length of data column and headers should match
        if len(hit_inform) != len(self._columns):
            raise ValueError(
                "The number of data of blastn hit differs from the number of headers specified. Headers {!r}, data {!r}".format(self._columns, hit_inform))

        # should NOT contain any sequence data column
        if any(x in self._columns for x in BLASTN_SEQ_COLUMNS):
            raise ValueError(
                "BlastnTSVHit class does not handle blastn output with sequence information (sseq/qseq specifier).")

    def __str__(self):
        return self._hit_str

    def __repr__(self):
        """
        Overwrite to have a better printing results with pprint (easier for debug).
        """
        return str(self)

    @property
    def columns(self):
        """
        Returns column names
        :return: _columns
        """
        return self._columns

    @property
    def has_data(blastn_data_specifier):
        """
        Check whether blastn hit object has specific data
        :param blastn_data_specifier: a specifier for blastn outfmt 6
        :return: True if has data, False otherwise
        """
        return blastn_data_specifier in self._columns
