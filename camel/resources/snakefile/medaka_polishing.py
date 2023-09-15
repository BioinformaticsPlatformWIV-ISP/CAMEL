from pathlib import Path

SNAKEFILE_MEDAKA_POLISHING = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_medaka = Path('medaka', '{assembly_type}')
INPUT_ASSEMBLY_FASTA = _dir_medaka / 'input' / 'fasta.io'
OUTPUT_ASSEMBLY_FASTA = _dir_medaka / 'stitch' / 'fasta.io'
OUTPUT_ASSEMBLY_MAPPING_INFORMS = _dir_medaka / 'minimap2' / 'informs.io'
OUTPUT_ASSEMBLY_MAPPING_RATE_INFORMS = _dir_medaka / 'samtools_flagstat' / 'informs.io'
OUTPUT_ASSEMBLY_DEPTH_INFORMS = _dir_medaka / 'samtools_depth' / 'informs.io'
OUTPUT_ASSEMBLY_REPORT = _dir_medaka / 'report' / 'html.io'
OUTPUT_ASSEMBLY_REPORT_EMPTY = _dir_medaka / 'report' / 'html_empty.io'
OUTPUT_ASSEMBLY_SUMMARY = _dir_medaka / 'summary' / 'summary_out.tsv'
