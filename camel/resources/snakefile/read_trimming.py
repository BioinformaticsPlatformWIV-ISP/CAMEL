import os


FOLDER_TRIMMING = 'read_trimming'
OUTPUT_READ_TRIMMING_REPORT = os.path.join('read_trimming', 'report', 'html.io')
OUTPUT_READ_TRIMMING_SUMMARY = os.path.join('read_trimming', 'summary', 'summary_out.tsv')
OUTPUT_READ_TRIMMING_READS_PE = os.path.join('read_trimming', 'trimmomatic', 'fastq-pe.io')
OUTPUT_READ_TRIMMING_READS_SE_FWD = os.path.join('read_trimming', 'trimmomatic', 'fastq-se-forward.io')
OUTPUT_READ_TRIMMING_READS_SE_REV = os.path.join('read_trimming', 'trimmomatic', 'fastq-se-reverse.io')
OUTPUT_READ_TRIMMING_INFORMS = os.path.join('read_trimming', 'trimmomatic', 'informs.io')
OUTPUT_READ_TRIMMING_FASTQC_POST = os.path.join('read_trimming', 'fastqc-post', 'html.io')
OUTPUT_READ_TRIMMING_FASTQC_PRE = os.path.join('read_trimming', 'fastqc-pre', 'html.io')