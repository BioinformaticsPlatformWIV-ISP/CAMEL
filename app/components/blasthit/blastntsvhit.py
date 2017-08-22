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
        return self._column

    # TODO remove following pipeline specific function

    # def get_hit_id(self):
    # backward compatibility
    #     return self.hit_id

    # def __check_validity(self, hit_inform):
    #     """
    #     Check the validity of the blast hit & its columns
    #     :param hit_inform: the blast hit data in list
    #     """
    #     if len(hit_inform) != len(self.columns):
    #         raise ValueError("Number of data column differs from the number of headers.")

    #     for key in REQUIRED_COLUMNS:
    #         if key not in self.columns:
    #             raise ValueError("Required blast hit data column {!r} is missing.".format(key))

    #     for col_name in self.columns:
    #         if col_name not in SUPPORTED_COLUMNS:
    #             raise ValueError("Not supported blastn tsv output data {!r}.".format(col_name))

    # @property
    # def hit_id(self):
    #     """
    #     Generate a unique hit id with sstart, send, and sstrand information
    #     :return: hit id as string
    #     """
    # Example id: NODE_10_length_1700_cov_25.444118:25-983(-)
    #     if self.strand == 'minus':
    #         return self.sseqid + ":" + str(self.sstart) + "-" + str(self.send) + "(-)"
    #     else:
    #         return self.sseqid + ":" + str(self.sstart) + "-" + str(self.send) + "(+)"

    # @staticmethod
    # def get_location_from_hit_id(hit_id):
    #     """
    #     Given hit_id generated using above function, return location information
    #     :param hit_id: blastn hit id generated using self.get_hit_id function
    #     :return: sstart
    #     :return: send
    #     :return: strand
    #     """
    # Example id: NODE_10_length_1700_cov_25.444118:25-983(-)
    #     location = hit_id.split(":")[-1]
    #     m = re.search(r"^(\d+)-(\d+)\(([+-])\)$", location)
    #     start, end, strand = m.groups()
    #     start = int(start)
    #     end = int(end)
    #     return start, end, strand
