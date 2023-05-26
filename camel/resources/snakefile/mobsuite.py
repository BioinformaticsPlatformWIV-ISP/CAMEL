from pathlib import Path

DIR_MOB_SUITE = Path('mob_suite')
SNAKEFILE_MOB_SUITE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

# Input
INPUT_MOBSUITE_FASTA = DIR_MOB_SUITE / 'input' / 'fasta.io'

# Report and summary
OUTPUT_MOB_SUITE_INFORMS = DIR_MOB_SUITE / 'informs.io'
OUTPUT_MOB_SUITE_REPORT = DIR_MOB_SUITE / 'html.io'
OUTPUT_MOB_SUITE_REPORT_EMPTY = DIR_MOB_SUITE / 'html-empty.io'
OUTPUT_MOB_SUITE_SUMMARY = DIR_MOB_SUITE / 'summary_mob_suite.tsv'
