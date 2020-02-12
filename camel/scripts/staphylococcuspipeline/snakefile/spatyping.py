from pathlib import Path

SNAKEFILE_SPATYPING = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_spatyping = Path('spa_typing')
OUTPUT_SPATYPING_REPORT = _dir_spatyping / 'report' / 'html.io'
OUTPUT_SPATYPING_REPORT_EMPTY = _dir_spatyping / 'report' / 'html-empty.io'
OUTPUT_SPATYPING_SUMMARY = _dir_spatyping / 'summary_spatyping.tsv'
