from pathlib import Path

SNAKEFILE_ANTIVIRALS = Path(__file__).parent / f'{Path(__file__).stem}.smk'
OUTPUT_REPORT = Path('antivirals', 'report', 'html.io')
OUTPUT_REPORT_EMPTY = Path('antivirals', 'report', 'html-empty.io')
OUTPUT_SUMMARY = Path('antivirals', 'summary.tsv')
