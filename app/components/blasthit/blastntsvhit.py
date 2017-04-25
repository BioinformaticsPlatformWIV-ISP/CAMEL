class BlastnTSVHit(object):

    """
    Class to handle customized blastn tabular output

    Supported blastn outfmt 6 columns:
    'qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore strand qcovs qcovhsp'
    """
    # Note that although one would expect that qcovs and qcovhsp is float, in fact, they are reported as integer (as I
    # never see decimal digits reported). Indirectly, another package also think they are integer (see below). A direct
    # proof will need to check the source code.
    #
    # Reference: http://scikit-bio.org/docs/0.4.1/generated/skbio.io.format.blast6.html
    DEFAULT_COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch",
                       "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    INT_COLUMNS = ['length', 'mismatch', 'gapopen', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'qcovs', 'qcovhsp']
    FLOAT_COLUMNS = ['pident', 'bitscore']

    # TODO MOVED TO PIPELINE SPECIFIC CODE
    SEQEXTRACTION_COLUMNS = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "gaps",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "strand", "qcovs", "qcovhsp"
    ]

    def __init__(self, hit, columns=BlastnTSVHit.DEFAULT_COLUMNS):
        """
        Initialize
        :param hit: the txt string describe hit
        :param columns: specified data column names
        """
        hit_inform = hit.split("\t")
        self._columns = columns

        for idx, col_name in enumerate(columns):

            try:
                if col_name in BlastnTSVHit.INT_COLUMNS:
                    setattr(self, col_name, int(hit_inform[idx]))
                elif col_name in BlastnTSVHit.FLOAT_COLUMNS:
                    setattr(self, col_name, float(hit_inform[idx]))
                else:
                    setattr(self, col_name, hit_inform[idx])
            except ValueError:
                logging.error("ValueError handling column {!r} at pos {} with value {!r}.".format(
                    col_name, idx, hit[idx]))
                raise

        self._hit_str = hit

    def __str__(self):
        return self._hit_str

    def __repr__(self):
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
