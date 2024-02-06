from pathlib import Path

SNAKEFILE_POLISH_ASSEMBLY_SHORT = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_polish = Path('polish', 'short_reads', '{assembly_type}')
INPUT_ASSEMBLY_FASTA = _dir_polish / 'input' / 'fasta.io'
OUTPUT_POLISHING_FASTA = _dir_polish / 'polca' / 'fasta.io'
OUTPUT_POLISHING_FASTA_INDEX_POLYPOLISH = _dir_polish / 'polypolish' / 'fasta-index.io'
OUTPUT_POLYPOLISH_INFORMS = _dir_polish / 'polypolish' / 'informs.io'
OUTPUT_POLYPOLISH_FASTA = _dir_polish / 'polypolish' / 'fasta.io'
OUTPUT_POLISHING_FASTA_INDEX_POLCA = _dir_polish / 'polca' / 'fasta-index.io'
OUTPUT_POLCA_INFORMS = _dir_polish / 'polca' / 'informs.io'
