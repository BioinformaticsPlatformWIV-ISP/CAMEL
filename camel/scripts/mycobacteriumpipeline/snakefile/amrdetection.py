from pathlib import Path

SNAKEFILE_AMR = f'{Path(__file__).parent / Path(__file__).stem}.smk'
OUTPUT_AMR_REPORT = Path('amr') / 'report' / 'html.io'
OUTPUT_AMR_REPORT_EMPTY = Path('amr') / 'report' / 'html-empty.io'
OUTPUT_AMR_SUMMARY = Path('amr') / 'summary' / 'summary_out.tsv'
OUTPUT_INFORMS_CSQ = Path('amr') / 'csq' / 'informs.io'
