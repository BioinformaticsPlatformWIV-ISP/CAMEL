from pathlib import Path

SNAKEFILE_READ_SIMULATION = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_simulation = Path('read_simulation')

OUTPUT_SIMULATION_FASTQ = _dir_simulation / 'art' / 'fastq.io'
OUTPUT_SIMULATION_INFORMS = _dir_simulation / 'art' / 'informs.io'
