from pathlib import Path

DIR_BACMET = Path('bacmet')
SNAKEFILE_BACMET = f'{Path(__file__).parent / Path(__file__).stem}.smk'

# Report and summary
OUTPUT_BACMET_INFORMS = DIR_BACMET / 'blastp' / 'informs.io'
OUTPUT_BACMET_REPORT = DIR_BACMET / 'report' / 'html.io'
OUTPUT_BACMET_REPORT_EMPTY = DIR_BACMET / 'report' / 'html-empty.io'
OUTPUT_BACMET_SUMMARY = DIR_BACMET / 'summary_bacmet.tsv'

# Prodigal
OUTPUT_PRODIGAL_REPORT = DIR_BACMET / 'prodigal' / 'report' / 'html.io'
OUTPUT_PRODIGAL_REPORT_EMPTY = DIR_BACMET / 'prodigal' / 'report' / 'html-empty.io'
