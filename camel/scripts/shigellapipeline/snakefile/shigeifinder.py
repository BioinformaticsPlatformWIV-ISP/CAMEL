from pathlib import Path


SNAKEFILE_SHIGEIFINDER = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_shigeifinder = Path('shigeifinder')
OUTPUT_SHIGEIFINDER_REPORT = _dir_shigeifinder / 'report' / 'html.io'
OUTPUT_SHIGEIFINDER_REPORT_EMPTY = _dir_shigeifinder / 'report' / 'html-empty.io'
OUTPUT_SHIGEIFINDER_SUMMARY = _dir_shigeifinder / 'summary_shigeifinder.tsv'
OUTPUT_SHIGEIFINDER_INFORMS = _dir_shigeifinder / 'informs.io'
