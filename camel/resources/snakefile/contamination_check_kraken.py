from pathlib import Path

SNAKEFILE_CONTAMINATION_CHECK_KRAKEN = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_out = Path('contamination_check')
OUTPUT_CONTAMINATION_CHECK_REPORT = _dir_out / 'report' / 'html.io'
OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY = _dir_out / 'report' / 'html-empty.io'
OUTPUT_CONTAMINATION_CHECK_INFORMS = _dir_out / 'kraken' / 'informs-contamination.io'
OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS = _dir_out / 'kraken' / 'informs.io'
OUTPUT_CONTAMINATION_SUMMARY = _dir_out / 'summary' / 'summary_out.tsv'
