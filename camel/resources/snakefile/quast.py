from pathlib import Path

SNAKEFILE_QUAST = f'{Path(__file__).parent / Path(__file__).stem}.smk'
OUTPUT_QUAST_REPORT = Path('quast') / 'report' / 'html.io'
OUTPUT_QUAST_SUMMARY = Path('quast') / 'report' / 'summary_quast.tsv'
OUTPUT_QUAST_INFORMS = Path('quast') / 'output' / 'informs.io'
