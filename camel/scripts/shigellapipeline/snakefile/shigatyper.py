from pathlib import Path


SNAKEFILE_SHIGATYPER = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_shigatyper = Path('shigatyper')
OUTPUT_SHIGATYPER_REPORT = _dir_shigatyper/ 'report' / 'html.io'
OUTPUT_SHIGATYPER_REPORT_EMPTY = _dir_shigatyper / 'report' / 'html-empty.io'
OUTPUT_SHIGATYPER_SUMMARY = _dir_shigatyper / 'summary_shigatyper.tsv'
