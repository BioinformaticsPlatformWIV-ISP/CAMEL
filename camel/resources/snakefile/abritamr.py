from pathlib import Path

SNAKEFILE_ABRITAMR = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_abritamr = Path('abritamr')
OUTPUT_MATCHES_ABRITAMR = _dir_abritamr / 'abritamr_output_matches.io'
OUTPUT_PARTIALS_ABRITAMR = _dir_abritamr / 'abritamr_output_partials.io'
OUTPUT_COMBINED_ABRITAMR = _dir_abritamr / 'abritamr_output_combined.io'
OUTPUT_QC_ABRITAMR = _dir_abritamr / 'qc_file.txt'
OUTPUT_ABRITAMR_RUN_INFORMS = _dir_abritamr / 'informs_run.io'
OUTPUT_REPORT_ABRITAMR = _dir_abritamr / 'abritamr_output_report.io'
OUTPUT_REPORT_ABRITAMR_INFORMS = _dir_abritamr / 'informs_report.io'
OUTPUT_ABRITAMR_REPORT = _dir_abritamr / 'html.io'
OUTPUT_ABRITAMR_REPORT_EMPTY = _dir_abritamr / 'html-empty.io'
OUTPUT_ABRITAMR_SUMMARY = _dir_abritamr / 'summary_out.tsv'
OUTPUT_ABRITAMR_SUMMARY_JSON = _dir_abritamr / 'summary_out.json'
