from pathlib import Path

DIR_RESFINDER4 = Path('resfinder4')
SNAKEFILE_RESFINDER4 = f'{Path(__file__).parent / Path(__file__).stem}.smk'

# Report and summary
OUTPUT_RESFINDER4_INFORMS = DIR_RESFINDER4 / 'informs.io'
OUTPUT_RESFINDER4_REPORT = DIR_RESFINDER4 / 'html.io'
OUTPUT_RESFINDER4_REPORT_EMPTY = DIR_RESFINDER4 / 'html-empty.io'
OUTPUT_RESFINDER4_SUMMARY = DIR_RESFINDER4 / 'summary_resfinder.tsv'
