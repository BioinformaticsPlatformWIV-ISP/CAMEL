from pathlib import Path

SNAKEFILE_SUBSPECIES = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_subspecies = Path('subspecies_identification')
OUTPUT_SPECIES_REPORT = _dir_subspecies / 'report' / 'html-species.io'
OUTPUT_SPECIES_REPORT_EMPTY = _dir_subspecies / 'report' / 'html-species-empty.io'
OUTPUT_SUBSPECIES_REPORT = _dir_subspecies / 'report' / 'html-subspecies.io'
OUTPUT_SUBSPECIES_REPORT_EMPTY = _dir_subspecies / 'report' / 'html-subspecies-empty.io'
OUTPUT_SUBSPECIES_SUMMARY = _dir_subspecies / 'summary_out.tsv'
