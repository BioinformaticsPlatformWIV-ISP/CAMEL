from pathlib import Path


SNAKEFILE_POINTFINDER = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_pointfinder = Path('pointfinder')
INPUT_POINTFINDER_FASTA = _dir_pointfinder / 'input' / 'fasta.io'
OUTPUT_POINTFINDER_REPORT = _dir_pointfinder / 'html.io'
OUTPUT_POINTFINDER_REPORT_EMPTY = _dir_pointfinder / 'html-empty.io'
OUTPUT_POINTFINDER_SUMMARY = _dir_pointfinder / 'summary_out.tsv'
OUTPUT_POINTFINDER_INFORMS = _dir_pointfinder / 'informs.io'
