import os

GENE_DETECTION_FASTA = os.path.join('gene_detection', '{db}', 'db_manager', 'fasta.io'),
GENE_DETECTION_FASTA_CLUSTERED = os.path.join('gene_detection', '{db}', 'db_manager', 'fasta-clust.io'),

INPUT_GENE_DETECTION_FASTA = os.path.join('gene_detection', '{db}', 'input', 'fasta.io')
INPUT_GENE_DETECTION_FASTQ_PE = os.path.join('gene_detection', '{db}', 'input', 'fastq-pe.io')

OUTPUT_GENE_DETECTION_HITS_BLAST = os.path.join('gene_detection', '{db}', 'alignment_extraction', 'blast-hits.io')
OUTPUT_GENE_DETECTION_HITS_SRST2 = os.path.join('gene_detection', '{db}', 'srst2', 'srst2-hits.io')

OUTPUT_GENE_DETECTION_TSV_BLAST = os.path.join('gene_detection', '{db}', 'hit_filtering', 'tsv-filtered.io')
OUTPUT_GENE_DETECTION_TSV_SRST2 = os.path.join('gene_detection', '{db}', 'hit_extraction', 'tsv-srst2.io')

OUTPUT_GENE_DETECTION_REPORT = os.path.join('gene_detection', '{db}', 'report', 'html.io')
OUTPUT_GENE_DETECTION_REPORT_EMPTY = os.path.join('gene_detection', '{db}', 'report', 'html-empty.io')
OUTPUT_GENE_DETECTION_SUMMARY = os.path.join('gene_detection', '{db}', 'report', 'summary_out.tsv')
