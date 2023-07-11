from pathlib import Path


SNAKEFILE_CONFINDR = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_confindr = Path('confindr')
OUTPUT_CONFINDR_INFORMS = _dir_confindr / 'informs.io'
OUTPUT_CONFINDR_REPORT = _dir_confindr / 'report' / 'html.io'
OUTPUT_CONFINDR_REPORT_EMPTY = _dir_confindr / 'report' / 'html-empty.io'
OUTPUT_CONFINDR_SUMMARY = _dir_confindr / 'summary' / 'summary_confindr.tsv'
