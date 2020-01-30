from pathlib import Path

SNAKEFILE_SEROTYPE = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_serotype = Path('serotype_detection')
OUTPUT_VAL_SEROTYPE = _dir_serotype / 'serotype_detection' / 'val-sero.io'
OUTPUT_SEROTYPE_REPORT = _dir_serotype / 'serotype_detection' / 'html.io'
OUTPUT_SEROTYPE_SUMMARY = _dir_serotype / 'serotype_detection' / 'summary_out.tsv'
