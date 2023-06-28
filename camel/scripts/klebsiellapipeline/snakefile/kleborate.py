from pathlib import Path

DIR_KLEBORATE = Path('kleborate')
SNAKEFILE_KLEBORATE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

# Report and summary
OUTPUT_KLEBORATE_INFORMS = DIR_KLEBORATE / 'informs.io'
OUTPUT_KLEBORATE_REPORT = DIR_KLEBORATE / 'html.io'
OUTPUT_KLEBORATE_REPORT_EMPTY = DIR_KLEBORATE / 'html-empty.io'
OUTPUT_KLEBORATE_SUMMARY = DIR_KLEBORATE / 'summary_kleborate.tsv'
