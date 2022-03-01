from pathlib import Path

SNAKEFILE_DOWNSAMPLING = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_downsampling = Path('downsampling')
OUTPUT_DOWNSAMPLING_REPORT = _dir_downsampling / 'html.io'
OUTPUT_DOWNSAMPLING_FASTQ_PE = _dir_downsampling / 'fastq-pe.io'
