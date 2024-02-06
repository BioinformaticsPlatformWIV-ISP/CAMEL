from pathlib import Path

SNAKEFILE_POLISHING = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_polishing = Path('polishing', '{assembly_type}')
INPUT_ASSEMBLY_FASTA = _dir_polishing / 'input' / 'fasta.io'
OUTPUT_POLISHING_FASTA = _dir_polishing / 'polca' / 'fasta.io'
OUTPUT_POLISHING_FASTA_INDEX_POLYPOLISH = _dir_polishing / 'polypolish' / 'fasta-index.io'
OUTPUT_POLYPOLISH_INFORMS = _dir_polishing / 'polypolish' / 'informs.io'
OUTPUT_POLYPOLISH_FASTA = _dir_polishing / 'polypolish' / 'fasta.io'
OUTPUT_POLISHING_FASTA_INDEX_POLCA = _dir_polishing / 'polca' / 'fasta-index.io'
OUTPUT_POLCA_INFORMS = _dir_polishing / 'polca' / 'informs.io'
