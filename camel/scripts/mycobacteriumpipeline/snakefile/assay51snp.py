from pathlib import Path

SNAKEFILE_51SNP = f'{Path(__file__).parent / Path(__file__).stem}.smk'
OUTPUT_51SNP_REPORT = Path('51snp') / 'html.io'
OUTPUT_51SNP_REPORT_EMPTY = Path('51snp') / 'html-empty.io'
OUTPUT_51SNP_SUMMARY = Path('51snp') / 'summary_out.tsv'
