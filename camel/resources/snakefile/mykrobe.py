from pathlib import Path

SNAKEFILE_MYKROBE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_mykrobe = Path('mykrobe')
OUTPUT_MYKROBE_REPORT = _dir_mykrobe / 'report' / 'html.io'
OUTPUT_MYKROBE_REPORT_EMPTY = _dir_mykrobe / 'report' / 'html-empty.io'
OUTPUT_MYKROBE_SUMMARY = _dir_mykrobe / 'summary_mykrobe.tsv'
OUTPUT_MYKROBE_INFORMS = _dir_mykrobe / 'informs.io'
