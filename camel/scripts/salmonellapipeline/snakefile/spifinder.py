from pathlib import Path

SNAKEFILE_SPIFINDER= f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_spifinder = Path('spifinder')
OUTPUT_SPIFINDER_FASTQ = _dir_spifinder / 'spifinder_fastq'  / 'spifinder_output.io'
OUTPUT_SPIFINDER_FASTA = _dir_spifinder / 'spifinder_fasta' / 'spifinder_output.io'
OUTPUT_SPIFINDER_REPORT = _dir_spifinder / 'html.io'
OUTPUT_SPIFINDER_FASTQ_INFORMS = _dir_spifinder / 'spifinder_fastq' / 'informs.io'
OUTPUT_SPIFINDER_FASTA_INFORMS = _dir_spifinder / 'spifinder_fasta' / 'informs.io'
OUTPUT_SPIFINDER_REPORT_EMPTY = _dir_spifinder / 'html-empty.io'
OUTPUT_SPIFINDER_SUMMARY = _dir_spifinder / 'summary_out.tsv'
OUTPUT_SPIFINDER_SUMMARY_JSON = _dir_spifinder / 'summary_out.json'
OUTPUT_SPIFINDER_DOC = _dir_spifinder / 'spifinder_function_category.tsv'
OUTPUT_SPIFINDER_SUMMARY_IO = _dir_spifinder / 'summary_out.io'
