from pathlib import Path

SNAKEFILE_SNPIT = f'{Path(__file__).parent / Path(__file__).stem}.smk'
OUTPUT_SNPIT_REPORT = Path('snpit') / 'html.io'
OUTPUT_SNPIT_REPORT_EMPTY = Path('snpit') / 'html-empty.io'
OUTPUT_SNPIT_SUMMARY = Path('snpit') / 'summary_out.tsv'
OUTPUT_SNPIT_INFORMS = Path('snpit') / 'informs.io'
