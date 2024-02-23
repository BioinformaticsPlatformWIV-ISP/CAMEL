from pathlib import Path

SNAKEFILE_SNP_LINEAGE = f'{Path(__file__).parent / Path(__file__).stem}.smk'
OUTPUT_SNP_LINEAGE_REPORT = Path('snp_lineage') / 'html.io'
OUTPUT_SNP_LINEAGE_REPORT_EMPTY = Path('snp_lineage') / 'html-empty.io'
OUTPUT_SNP_LINEAGE_SUMMARY = Path('snp_lineage') / 'summary_out.tsv'
OUTPUT_SNP_LINEAGE_INFORMS = Path('snp_lineage') / 'informs.io'
