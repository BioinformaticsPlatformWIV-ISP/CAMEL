from pathlib import Path


SNAKEFILE_ASSEMBLY_CANU = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_canu = Path('assembly_canu')
OUTPUT_ASSEMBLY_FASTA = _dir_canu / 'filtering' / 'fasta.io'
OUTPUT_ASSEMBLY_INFORMS = _dir_canu / 'canu' / 'informs.io'
OUTPUT_ASSEMBLY_FILTERING_INFORMS = _dir_canu / 'filtering' / 'informs.io'
OUTPUT_ASSEMBLY_REPORT = _dir_canu / 'report' / 'html.io'
OUTPUT_ASSEMBLY_SUMMARY = _dir_canu / 'summary' / 'summary_out.tsv'
