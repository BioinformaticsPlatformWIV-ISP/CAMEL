from pathlib import Path

SNAKEFILE_LREFINDER = f'{Path(__file__).parent / Path(__file__).stem}.smk'
OUTPUT_LREFINDER_REPORT = Path('lrefinder') / 'html.io'
OUTPUT_LREFINDER_REPORT_EMPTY = Path('lrefinder') / 'html_empty.io'
OUTPUT_LREFINDER_INFORMS = Path('lrefinder') / 'informs.io'
OUTPUT_LREFINDER_SUMMARY = Path('lrefinder') / 'summary.tsv'
