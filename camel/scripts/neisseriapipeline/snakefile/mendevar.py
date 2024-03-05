from pathlib import Path


SNAKEFILE_MENDEVAR = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_gmats = Path('mendevar')
OUTPUT_MENDEVAR_REPORT = _dir_gmats / 'report' / 'html.io'
OUTPUT_MENDEVAR_REPORT_EMPTY = _dir_gmats / 'report' / 'html-empty.io'
OUTPUT_MENDEVAR_SUMMARY = _dir_gmats / 'summary' / 'summary_mendevar.tsv'
