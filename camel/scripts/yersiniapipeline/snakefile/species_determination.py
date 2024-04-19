from pathlib import Path

SNAKEFILE_SPECIES_DETERMINATION = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_species = Path('species_determination')
OUTPUT_SPECIES_DETERMINATION_REPORT = _dir_species / 'html.io'
OUTPUT_SPECIES_DETERMINATION_REPORT_EMPTY = _dir_species / 'html-empty.io'
OUTPUT_SPECIES_DETERMINATION_INFORMS = _dir_species / 'informs.io'
OUTPUT_SPECIES_DETERMINATION_SUMMARY = _dir_species / 'summary_species_determination.tsv'
OUTPUT_SPECIES_DETERMINATION_TSV = _dir_species / 'tsv_species_determination.io'
