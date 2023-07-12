from pathlib import Path

SNAKEFILE_SEROGROUP_DETERMINATION = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_serogroup = Path('serogroup_determination')
OUTPUT_SEROGROUP_DETERMINATION_LEGACY_REPORT = _dir_serogroup / 'legacy' / 'html.io'
OUTPUT_SEROGROUP_DETERMINATION_LEGACY_REPORT_EMPTY = _dir_serogroup / 'legacy' / 'html-empty.io'
OUTPUT_SEROGROUP_DETERMINATION_REPORT = _dir_serogroup / 'capsule' / 'html.io'
OUTPUT_SEROGROUP_DETERMINATION_INFORMS = _dir_serogroup / 'capsule' / 'informs.io'
OUTPUT_SEROGROUP_DETERMINATION_REPORT_EMPTY = _dir_serogroup / 'capsule' / 'html-empty.io'
OUTPUT_SEROGROUP_DETERMINATION_SUMMARY = _dir_serogroup / 'summary_out.tsv'
