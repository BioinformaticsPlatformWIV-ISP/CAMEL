from pathlib import Path

SNAKEFILE_ASSEMBLY_SPADES = f'{Path(__file__).parent / Path(__file__).stem}.smk'
OUTPUT_ASSEMBLY_FASTA = Path('assembly', 'spades', 'fasta.io')
OUTPUT_ASSEMBLY_INFORMS = Path('assembly', 'spades', 'informs.io')
