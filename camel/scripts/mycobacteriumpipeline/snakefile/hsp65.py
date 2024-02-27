from pathlib import Path

SNAKEFILE_HSP65 = f'{Path(__file__).parent / Path(__file__).stem}.smk'
OUTPUT_HSP65_REPORT = Path('hsp65') / 'html.io'
OUTPUT_HSP65_REPORT_EMPTY = Path('hsp65') / 'html-empty.io'
OUTPUT_HSP65_SUMMARY = Path('hsp65') / 'summary.tsv'
