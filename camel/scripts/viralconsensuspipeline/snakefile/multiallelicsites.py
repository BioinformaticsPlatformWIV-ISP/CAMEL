from pathlib import Path

SNAKEFILE_MULTI_ALLELIC = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_multi_allelic = Path('multi_allelic')
OUTPUT_MULTI_ALLELIC_REPORT = _dir_multi_allelic / 'report' / 'html.io'
OUTPUT_MULTI_ALLELIC_SUMMARY = _dir_multi_allelic / 'report' / 'summary.tsv'
