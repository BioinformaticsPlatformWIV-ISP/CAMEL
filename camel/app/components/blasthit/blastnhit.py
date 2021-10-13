import re
from typing import Union, Tuple


class BlastnHit(object):

    def __init__(self, qseqid: Union[str, int, float] = None, qgi: Union[str, int, float] = None, qacc: Union[str, int, float] = None,
                 qaccver: Union[str, int, float] = None, qlen: Union[str, int, float] = None, sseqid: Union[str, int, float] = None,
                 sallseqid: Union[str, int, float] = None, sgi: Union[str, int, float] = None, sallgi: Union[str, int, float] = None,
                 sacc: Union[str, int, float] = None, saccver: Union[str, int, float] = None, sallacc: Union[str, int, float] = None,
                 slen: Union[str, int, float] = None, qstart: Union[str, int, float] = None, qend: Union[str, int, float] = None,
                 sstart: Union[str, int, float] = None, send: Union[str, int, float] = None, qseq: Union[str, int, float] = None,
                 sseq: Union[str, int, float] = None, evalue: Union[str, int, float] = None, bitscore: Union[str, int, float] = None,
                 score: Union[str, int, float] = None, length: Union[str, int, float] = None, pident: Union[str, int, float] = None,
                 nident: Union[str, int, float] = None, mismatch: Union[str, int, float] = None, positive: Union[str, int, float] = None,
                 gapopen: Union[str, int, float] = None, gaps: Union[str, int, float] = None, ppos: Union[str, int, float] = None,
                 frames: Union[str, int, float] = None, qframe: Union[str, int, float] = None, sframe: Union[str, int, float] = None,
                 btop: Union[str, int, float] = None, staxids: Union[str, int, float] = None, sscinames: Union[str, int, float] = None,
                 scomnames: Union[str, int, float] = None, sblastnames: Union[str, int, float] = None, sskingdoms: Union[str, int, float] = None,
                 stitle: Union[str, int, float] = None, sstrand: Union[str, int, float] = None, salltitles: Union[str, int, float] = None,
                 qcovs: Union[str, int, float] = None, qcovhsp: Union[str, int, float] = None):
        """
        Initializes the hit object
        :param qseqid: Query Seq-id
        :param qgi: Query GI
        :param qacc: Query accesion
        :param qaccver: Query accesion.version
        :param qlen: Query sequence length
        :param sseqid: Subject Seq-id
        :param sallseqid: All subject Seq-id(s), separated by a ‘;’
        :param sgi: Subject GI
        :param sallgi: All subject GIs
        :param sacc: Subject accesion
        :param saccver: Subject accesion.version
        :param sallacc: All subject accesions
        :param slen: Subject sequence length
        :param qstart: Start of alignment in query
        :param qend: End of alignment in query
        :param sstart: Start of alignment in subject
        :param send: End of alignment in subject
        :param qseq: Aligned part of query sequence
        :param sseq: Aligned part of subject sequence
        :param evalue: Expect value
        :param bitscore: Bit score
        :param score: Raw score
        :param length: Alignment length
        :param pident: Percent of identical matches
        :param nident: Number of identical matches
        :param mismatch: Number of mismatches
        :param positive: Number of positive-scoring matches
        :param gapopen: Number of gap openings
        :param gaps: Total number of gaps
        :param ppos: Percentage of positive-scoring matches
        :param frames: Query and subject frames separated by a ‘/’
        :param qframe: Query frame
        :param sframe: Subject frame
        :param btop: Blast traceback operations (BTOP)
        :param staxids: Unique Subject Taxonomy ID(s), separated by a ‘;’ (in numerical order).
        :param sscinames: Unique Subject Scientific Name(s), separated by a ‘;’
        :param scomnames: Unique Subject Common Name(s), separated by a ‘;’
        :param sblastnames: unique Subject Blast Name(s), separated by a ‘;’ (in alphabetical order)
        :param sskingdoms: unique Subject Super Kingdom(s), separated by a ‘;’ (in alphabetical order)
        :param stitle: Subject Title
        :param sstrand: Subject Strand
        :param salltitles: All Subject Title(s), separated by a ‘<>’
        :param qcovs: Query Coverage Per Subject
        :param qcovhsp: Query Coverage Per HSP
        """
        self.qseqid = str(qseqid) if qseqid else None
        self.qgi = int(qgi) if qgi else None
        self.qacc = str(qacc) if qacc else None
        self.qaccver = str(qaccver) if qaccver else None
        self.qlen = int(qlen) if qlen else None
        self.sseqid = str(sseqid) if sseqid else None
        self.sallseqid = str(sallseqid) if sallseqid else None
        self.sgi = int(sgi) if sgi else None
        self.sallgi = int(sallgi) if sallgi else None
        self.sacc = str(sacc) if sacc else None
        self.saccver = str(saccver) if saccver else None
        self.sallacc = str(sallacc) if sallacc else None
        self.slen = int(slen) if slen else None
        self.qstart = int(qstart) if qstart else None
        self.qend = int(qend) if qend else None
        self.sstart = int(sstart) if sstart else None
        self.send = int(send) if send else None
        self.qseq = str(qseq) if qseq else None
        self.sseq = str(sseq) if sseq else None
        self.evalue = float(evalue) if evalue else None
        self.bitscore = float(bitscore) if bitscore else None
        self.score = int(score) if score else None
        self.length = int(length) if length else None
        self.pident = float(pident) if pident else None
        self.nident = int(nident) if nident else None
        self.mismatch = int(mismatch) if mismatch else None
        self.positive = int(positive) if positive else None
        self.gapopen = int(gapopen) if gapopen else None
        self.gaps = int(gaps) if gaps else None
        self.ppos = float(ppos) if ppos else None
        self.frames = str(frames) if frames else None
        self.qframe = int(qframe) if qframe else None
        self.sframe = int(sframe) if sframe else None
        self.btop = str(btop) if btop else None
        self.staxids = str(staxids) if staxids else None
        self.sscinames = str(sscinames) if sscinames else None
        self.scomnames = str(scomnames) if scomnames else None
        self.sblastnames = str(sblastnames) if sblastnames else None
        self.sskingdoms = str(sskingdoms) if sskingdoms else None
        self.stitle = str(stitle) if stitle else None
        self.sstrand = str(sstrand) if sstrand else None
        self.salltitles = str(salltitles) if salltitles else None
        self.qcovs = int(qcovs) if qcovs else None
        self.qcovhsp = int(qcovhsp) if qcovhsp else None

    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        return '\n'.join([f'{p}: {getattr(self, p)}' for p in self.__dict__.keys() if getattr(self, p) is not None and p not in {'sseq', 'qseq'}])

    def __repr__(self):
        return '\n'.join([f'{p}: {getattr(self, p)}' for p in self.__dict__.keys() if getattr(self, p) is not None and p not in {'sseq', 'qseq'}])

    def __eq__(self, other):
        if isinstance(other, BlastnHit):
            return other.get_hit_id() == self.get_hit_id()
        else:
            return False

    def get_hit_id(self) -> str:
        """
        Generate a unique hit id with sstart, send, and sstrand information
        :return: hit id as string
        """
        # Example id: NODE_10_length_1700_cov_25.444118:25-983(-)
        if self.sstrand == 'minus':
            return f'{self.sseqid}:{self.sstart}-{self.send}(-)'
        else:
            return f'{self.sseqid}:{self.sstart}-{self.send}(+)'

    @staticmethod
    def get_location_from_hit_id(hit_id: str) -> Tuple[int, int, str]:
        """
        Given the hit_id generated using the get_hit_id function, return the location information
        :param hit_id: blastn hit id generated using self.get_hit_id function
        :return: sstart
        :return: send
        :return: sstrand
        """
        # Example id: NODE_10_length_1700_cov_25.444118:25-983(-)
        location = hit_id.split(":")[-1]
        m = re.search(r"^(\d+)-(\d+)\(([+-])\)$", location)
        start, end, sstrand = m.groups()
        start = int(start)
        end = int(end)
        return start, end, sstrand

    @staticmethod
    def get_location_from_bedtool_seqid(seqid: str) -> Tuple[int, int, str]:
        """
        Given seqid generated by Bedtools getfasta function, return location information
        :param seqid: sequence id generated by Bedtools getfasta function
        :return: sstart
        :return: send
        :return: sstrand
        """
        # Example id: NODE_10_length_1700_cov_25.444118:25-983(-)
        #
        # NOTE that to extract the exact section of matching sequence, start position feed to Bedtools getfasta function
        # has been modified by 1 base (compare to the original blastn hit)
        location = seqid.split(":")[-1]
        m = re.search(r"^(\d+)-(\d+)\(([+-])\)$", location)
        start, end, sstrand = m.groups()
        start = int(start) + 1
        end = int(end)
        return start, end, sstrand
