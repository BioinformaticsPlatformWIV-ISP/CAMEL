from pathlib import Path

SNAKEFILE_POLISHING = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_polishing = Path('polishing', '{assembly_type}')
INPUT_ASSEMBLY_FASTA = _dir_polishing / 'input' / 'fasta.io'
INPUT_READS_FASTQ = Path('polishing') / 'input' / 'fastq.io'
OUTPUT_POLISHING_FASTA = _dir_polishing / 'fasta.io'
OUTPUT_POLISHING_REPORT = _dir_polishing / 'report' / 'html.io'
OUTPUT_POLISHING_REPORT_EMPTY = _dir_polishing / 'report' / 'html_empty.io'
OUTPUT_POLISHING_SUMMARY = _dir_polishing / 'summary' / 'summary_out.tsv'
