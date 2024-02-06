from pathlib import Path

SNAKEFILE_MEDAKA_POLISHING = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_medaka = Path('medaka', '{assembly_type}')
INPUT_ASSEMBLY_FASTA = _dir_medaka / 'input' / 'fasta.io'
OUTPUT_ASSEMBLY_FASTA = _dir_medaka / 'stitch' / 'fasta.io'
OUTPUT_POLISH_MEDAKA_INFORMS = _dir_medaka / 'stitch' / 'commands-stitch.io'
OUTPUT_ASSEMBLY_REPORT_EMPTY = _dir_medaka / 'report' / 'html_empty.io'
