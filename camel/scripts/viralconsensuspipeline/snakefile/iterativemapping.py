from pathlib import Path

SNAKEFILE_ITERATIVE_MAPPING = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_iterative_mapping = Path('iterative_mapping')
INPUT_FASTA_REF = _dir_iterative_mapping / 'input' / 'fasta.io'
INPUT_FASTQ = _dir_iterative_mapping / 'input' / 'fastq.io'
OUTPUT_ITERATIVE_MAPPING_REPORT = _dir_iterative_mapping / 'report' / 'html.io'
OUTPUT_ITERATIVE_MAPPING_SUMMARY = _dir_iterative_mapping / 'report' / 'summary_iterative_mapping.tsv'
OUTPUT_ITERATIVE_MAPPING_INFORMS = _dir_iterative_mapping / 'report' / 'informs.io'
OUTPUT_ITERATIVE_MAPPING_FASTA_CONSENSUS_FINAL = _dir_iterative_mapping / 'output' / 'fasta.io'
