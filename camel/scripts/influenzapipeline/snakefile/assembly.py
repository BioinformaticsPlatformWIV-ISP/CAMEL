from pathlib import Path

SNAKEFILE_ASSEMBLY = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_assembly = Path('assembly')
OUTPUT_ASSEMBLY_FASTA = _dir_assembly / 'fasta.io'
OUTPUT_ASSEMBLY_FASTA_BLAST_DB = _dir_assembly / 'blast_db.io'
