from pathlib import Path

SNAKEFILE_CSB_RD = f'{Path(__file__).parent / Path(__file__).stem}.smk'
OUTPUT_CSB_RD_REPORT = Path('csb_rd') / 'html.io'
OUTPUT_CSB_RD_REPORT_EMPTY = Path('csb_rd') / 'html-empty.io'
OUTPUT_CSB_RD_SUMMARY = Path('csb_rd') / 'summary.tsv'
