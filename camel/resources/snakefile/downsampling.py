from pathlib import Path

SNAKEFILE_DOWNSAMPLING = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_downsampling = Path('downsampling')
OUTPUT_DOWNSAMPLING_REPORT = _dir_downsampling / 'html.io'
OUTPUT_DOWNSAMPLING_SUMMARY = _dir_downsampling / 'summary_downsampling.tsv'
OUTPUT_DOWNSAMPLING_FASTQ = _dir_downsampling / 'fastq.io'
OUTPUT_DOWNSAMPLING_INFORMS = _dir_downsampling / 'downsampling' / 'informs_seqtk.io'
