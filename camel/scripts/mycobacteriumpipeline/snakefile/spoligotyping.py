from pathlib import Path

SNAKEFILE_SPOLIGOTYPING = f'{Path(__file__).parent / Path(__file__).stem}.smk'
OUTPUT_SPOLIGOTYPING_REPORT = Path('spoligotyping') / 'html.io'
OUTPUT_SPOLIGOTYPING_REPORT_EMPTY = Path('spoligotyping') / 'html-empty.io'
OUTPUT_SPOLIGOTYPING_SUMMARY = Path('spoligotyping') / 'summary_out.tsv'
OUTPUT_SPOLIGOTYPING_INFORMS = Path('spoligotyping') / 'informs.io'
