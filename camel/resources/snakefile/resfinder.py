from pathlib import Path

SNAKEFILE_RESFINDER = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_resfinder = Path('resfinder')
OUTPUT_RESFINDER_PHENO = _dir_resfinder / 'val-resfinder-pheno.io'
OUTPUT_RESFINDER_GENE = _dir_resfinder / 'val-resfinder-gene.io'
OUTPUT_RESFINDER_INFORMS = _dir_resfinder / 'informs.io'
OUTPUT_RESFINDER_REPORT = _dir_resfinder / 'html.io'
OUTPUT_RESFINDER_REPORT_EMPTY = _dir_resfinder / 'html-empty.io'
OUTPUT_RESFINDER_SUMMARY = _dir_resfinder / 'summary_out.tsv'
INPUT_RESFINDER_FASTA = _dir_resfinder / 'fasta.io'
