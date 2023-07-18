from pathlib import Path

SNAKEFILE_ASSEMBLY_FLYE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_flye = Path('assembly_flye')
OUTPUT_ASSEMBLY_FASTA = _dir_flye / 'medaka' / 'fasta.io'
OUTPUT_ASSEMBLY_INFORMS = _dir_flye / 'flye' / 'informs.io'
OUTPUT_ASSEMBLY_FILTERING_INFORMS = _dir_flye / 'filtering' / 'informs.io'
OUTPUT_ASSEMBLY_REPORT = _dir_flye / 'report' / 'html.io'
OUTPUT_ASSEMBLY_SUMMARY = _dir_flye / 'summary' / 'summary_out.tsv'
OUTPUT_ASSEMBLY_MAPPING_RATE_INFORMS = _dir_flye / 'samtools_flagstat' / 'informs.io'
OUTPUT_ASSEMBLY_DEPTH_INFORMS = _dir_flye / 'samtools_depth' / 'informs.io'
OUTPUT_ASSEMBLY_MAPPING_INFORMS = _dir_flye / 'minimap2' / 'informs.io'
