from pathlib import Path

SNAKEFILE_POLISHING = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_polishing = Path('polishing', '{assembly_type}')
INPUT_ASSEMBLY_FASTA = _dir_polishing / 'input' / 'fasta.io'
OUTPUT_POLISHING_FASTA = _dir_polishing / 'output' / 'fasta.io'
OUTPUT_POLISHING_REPORT = _dir_polishing / 'report' / 'html.io'
OUTPUT_POLISHING_REPORT_EMPTY = _dir_polishing / 'report' / 'html_empty.io'
OUTPUT_POLISHING_SUMMARY = _dir_polishing / 'summary' / 'summary_out.tsv'
OUTPUT_ASSEMBLY_FILTERING_INFORMS = _dir_polishing / 'filtering' / 'informs.io'
OUTPUT_ASSEMBLY_REPORT = _dir_polishing / 'report' / 'html.io'
OUTPUT_ASSEMBLY_SUMMARY = _dir_polishing / 'summary' / 'summary_out.tsv'
OUTPUT_ASSEMBLY_INFORMS = _dir_polishing / 'polca' / 'informs.io'
OUTPUT_ASSEMBLY_MAPPING_INFORMS = _dir_polishing / 'bowtie2' / 'informs.io'
OUTPUT_ASSEMBLY_DEPTH_INFORMS = _dir_polishing / 'samtools_depth' / 'informs.io'
OUTPUT_ASSEMBLY_NANOPORE_DEPTH_INFORMS = _dir_polishing / 'samtools_depth_nanopore' / 'informs.io'
OUTPUT_ASSEMBLY_NANOPORE_MAPPING_INFORMS = _dir_polishing / 'minimap2' / 'informs.io'
OUTPUT_ASSEMBLY_NANOPORE_MAPPING_RATE_INFORMS = _dir_polishing / 'samtools_flagstat' / 'informs.io'
