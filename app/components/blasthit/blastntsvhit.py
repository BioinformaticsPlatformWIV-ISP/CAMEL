import logging
import re

from app.components.blasthit import SUPPORT_COLUMNS, DEFAULT_COLUMNS, INT_COLUMNS, FLOAT_COLUMNS


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

    # TODO remove get_hit_id

    def __init__(self, hit, columns=DEFAULT_COLUMNS):
        """
        Initialize
        :param hit: the txt string describe hit
        """
        hit_inform = hit.split("\t")

        if len(hit_inform) != len(columns):
            raise ValueError("Number of data column differs from the number of headers.")

        for idx, col_name in enumerate(columns):

            if col_name not in SUPPORT_COLUMNS:
                raise ValueError("Supported blastn tsv output data {!r}.".format(col_name))

            try:
                if col_name in INT_COLUMNS:
                    setattr(self, col_name, int(hit_inform[idx]))
                elif col_name in FLOAT_COLUMNS:
                    setattr(self, col_name, float(hit_inform[idx]))
                else:
                    setattr(self, col_name, hit_inform[idx])
            except ValueError:
                logging.error("ValueError handling column {!r} at pos {} with value {!r}.".format(
                    col_name, idx, hit[idx]))
                raise

        self._hit_str = hit
        self._columns = columns

    def __str__(self):
        return self._hit_str

    def __repr__(self):
        return str(self)

    def get_hit_id(self):
        # backward compatibility
        return self.hit_id

    @property
    def columns(self):
        """
        Returns column names
        :return: _columns
        """
        return self._columns

    @property
    def hit_id(self):
        """
        Generate a uniq hit id with sstart, send, and sstrand information
        :return: hit id as string
        """
        # Example id: NODE_10_length_1700_cov_25.444118:25-983(-)
        if self.strand == 'minus':
            return self.sseqid + ":" + str(self.sstart) + "-" + str(self.send) + "(-)"
        else:
            return self.sseqid + ":" + str(self.sstart) + "-" + str(self.send) + "(+)"

    @staticmethod
    def get_location_from_hit_id(hit_id):
        """
        Given hit_id generated using above function, return location information
        :param hit_id: blastn hit id generated using self.get_hit_id function
        :return: sstart
        :return: send
        :return: strand
        """
        # Example id: NODE_10_length_1700_cov_25.444118:25-983(-)
        location = hit_id.split(":")[-1]
        m = re.search(r"^(\d+)-(\d+)\(([+-])\)$", location)
        start, end, strand = m.groups()
        start = int(start)
        end = int(end)
        return start, end, strand
