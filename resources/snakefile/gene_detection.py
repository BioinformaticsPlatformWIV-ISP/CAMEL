import os

INPUT_GENE_DETECTION_FASTA = os.path.join('gene_detection', '{db}', 'input', 'fasta.io')
INPUT_GENE_DETECTION_FASTQ_PE = os.path.join('gene_detection', '{db}', 'input', 'fastq-pe.io')

OUTPUT_GENE_DETECTION_REPORT = os.path.join('gene_detection', '{db}', 'report', 'html.io')
OUTPUT_GENE_DETECTION_REPORT_EMPTY = os.path.join('gene_detection', '{db}', 'report', 'html-empty.io')
OUTPUT_GENE_DETECTION_SUMMARY = os.path.join('gene_detection', '{db}', 'report', 'summary_out.tsv')
