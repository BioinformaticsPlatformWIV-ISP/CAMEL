from pathlib import Path

SNAKEFILE_SCRUBBING = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_scrubbing = Path('human_read_scrubbing')
INPUT_SCRUBBING_FASTQ = _dir_scrubbing / 'input' / 'fastq.io'
INPUT_SCRUBBING_FASTA = _dir_scrubbing / 'input' / 'fasta.io'
OUTPUT_SCRUBBING_INFORMS = _dir_scrubbing / 'scrubbing' / 'informs.io'
OUTPUT_SCRUBBING_REPORT = _dir_scrubbing / 'output' / 'html.io'
OUTPUT_SCRUBBING_SUMMARY = _dir_scrubbing / 'output' / 'summary_out.tsv'
OUTPUT_SCRUBBING_SUMMARY_JSON = _dir_scrubbing / 'output' / 'summary_out.json'
OUTPUT_SCRUBBING_FASTQ = _dir_scrubbing / 'output' / 'fastq.io'
OUTPUT_SCRUBBING_FASTA = _dir_scrubbing / 'output' / 'fasta.io'
