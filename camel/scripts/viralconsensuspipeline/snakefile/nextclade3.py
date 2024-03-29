from pathlib import Path

SNAKEFILE_NEXTCLADE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_nextclade = Path('nextclade')
OUTPUT_NEXTCLADE_SUBTYPE_REPORT = _dir_nextclade / 'subtype_determination' / 'report' / 'html.io'
OUTPUT_NEXTCLADE_SUBTYPE_REPORT_EMPTY = _dir_nextclade / 'subtype_determination' / 'report' / 'html-empty.io'
OUTPUT_NEXTCLADE_REPORT = _dir_nextclade / 'html.io'
OUTPUT_NEXTCLADE_REPORT_EMPTY = _dir_nextclade / 'html-empty.io'
OUTPUT_NEXTCLADE_SUMMARY = _dir_nextclade / 'summary_nextclade.tsv'
OUTPUT_NEXTCLADE_INFORMS = _dir_nextclade / 'informs.io'
