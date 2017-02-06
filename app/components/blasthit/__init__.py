SUPPORT_COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "gaps", "qstart", "qend",
                   "sstart", "send", "evalue", "bitscore", "strand", "qcovs", "qcovhsp", "sseq", "qseq"]


DEFAULT_COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch",
                   "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

SEQEXTRACTION_COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "gaps",
                         "qstart", "qend", "sstart", "send", "evalue", "bitscore", "strand", "qcovs", "qcovhsp"]

SEQEXTRACTION_COLUMNS_WITH_SEQS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "gaps", "qstart", "qend",
                                   "sstart", "send", "qseq", "sseq", "evalue", "bitscore", "strand", "qcovs", "qcovhsp"]

INT_COLUMNS = ['length', 'mismatch', 'gapopen', 'gaps', 'qstart', 'qend', 'sstart', 'send', 'qcovs', 'qcovhsp']

FLOAT_COLUMNS = ['pident', 'bitscore']

SEQ_COLUMNS = ['sseq', 'qseq']
