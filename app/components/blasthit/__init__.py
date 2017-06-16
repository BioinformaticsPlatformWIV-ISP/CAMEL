# Blastn output supported format specifiers grouped by data types for parsing and automatic converting blastn output
#
#  - 'qcovs' & 'qcovhsp': although one would expect that qcovs and qcovhsp is float, in fact, they are reported as integer (as I
#      never see decimal digits reported). Indirectly, another package also think they are integer (see below). A direct
#      proof will need to check the source code.
#
#      Reference: http://scikit-bio.org/docs/0.4.1/generated/skbio.io.format.blast6.html

BLASTN_INT_COLUMNS = ('qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'score', 'length', 'nident', 'mismatch',
                      'positive', 'gapopen', 'gaps', 'qcovs', 'qcovhsp')

BLASTN_FLOAT_COLUMNS = ('pident', 'ppos', 'bitscore')

BLASTN_SEQ_COLUMNS = ('sseq', 'qseq')
