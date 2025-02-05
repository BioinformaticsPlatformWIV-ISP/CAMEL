from pathlib import Path

SNAKEFILE_POLISH_ASSEMBLY_LONG = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_polish = Path('polish', 'long_reads', '{assembly_type}')
INPUT_ASSEMBLY_FASTA = _dir_polish / 'input' / 'fasta.io'
OUTPUT_ASSEMBLY_FASTA = _dir_polish / 'sequence' / 'fasta.io'
OUTPUT_POLISH_MEDAKA_INFORMS = _dir_polish / 'sequence' / 'commands-sequence.io'
OUTPUT_ASSEMBLY_REPORT_EMPTY = _dir_polish / 'report' / 'html_empty.io'
