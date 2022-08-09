from pathlib import Path

SNAKEFILE_BTYPER = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_btyper = Path('btyper')
OUTPUT_VAL_BTYPER = _dir_btyper / 'val-btyper.io'
OUTPUT_INFORMS_BTYPER = _dir_btyper / 'informs.io'
OUTPUT_BTYPER_REPORT = _dir_btyper / 'html.io'
OUTPUT_BTYPER_REPORT_EMPTY = _dir_btyper / 'html-empty.io'
OUTPUT_BTYPER_SUMMARY = _dir_btyper / 'summary_out.tsv'
