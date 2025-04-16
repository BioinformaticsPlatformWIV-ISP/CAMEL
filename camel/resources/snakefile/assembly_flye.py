from pathlib import Path

SNAKEFILE_ASSEMBLY_FLYE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

OUTPUT_ASSEMBLY_FASTA = Path('assembly', 'flye', 'fasta.io')
OUTPUT_ASSEMBLY_INFORMS = Path('assembly', 'flye', 'informs.io')
