from pathlib import Path


SNAKEFILE_GMATS = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_gmats = Path('gmats')
OUTPUT_GMATS_REPORT = _dir_gmats / 'report' / 'html.io'
OUTPUT_GMATS_REPORT_EMPTY = _dir_gmats / 'report' / 'html-empty.io'
OUTPUT_GMATS_SUMMARY = _dir_gmats / 'summary' / 'summary_gmats.tsv'
