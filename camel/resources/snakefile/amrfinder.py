from pathlib import Path

DIR_AMRFINDER = Path('amrfinder')
SNAKEFILE_AMRFINDER = f'{Path(__file__).parent / Path(__file__).stem}.smk'

# Input
INPUT_AMRFINDER_FASTA = DIR_AMRFINDER / 'fasta.io'

# Report and summary
OUTPUT_AMRFINDER_INFORMS = DIR_AMRFINDER / 'informs.io'
OUTPUT_AMRFINDER_REPORT = DIR_AMRFINDER / 'html.io'
OUTPUT_AMRFINDER_REPORT_EMPTY = DIR_AMRFINDER / 'html-empty.io'
OUTPUT_AMRFINDER_SUMMARY = DIR_AMRFINDER / 'summary_amrfinder.tsv'
