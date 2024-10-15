from pathlib import Path
from typing import Dict, Any, List

SNAKEFILE_ABRITAMR = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_abritamr = Path('abritamr')
OUTPUT_MATCHES_ABRITAMR = _dir_abritamr / 'abritamr_output_matches.io'
OUTPUT_PARTIALS_ABRITAMR = _dir_abritamr / 'abritamr_output_partials.io'
OUTPUT_AMRFINDER_ABRITAMR = _dir_abritamr / 'abritamr_output_amrfinder.io'
OUTPUT_COMBINED_ABRITAMR = _dir_abritamr / 'abritamr_output_combined.io'
OUTPUT_QC_ABRITAMR = _dir_abritamr / 'qc_file.txt'
OUTPUT_ABRITAMR_RUN_INFORMS = _dir_abritamr / 'informs_run.io'
OUTPUT_REPORT_ABRITAMR = _dir_abritamr / 'abritamr_output_report.io'
OUTPUT_REPORT_ABRITAMR_INFORMS = _dir_abritamr / 'informs_report.io'
OUTPUT_ABRITAMR_REPORT = _dir_abritamr / 'html.io'
OUTPUT_ABRITAMR_REPORT_EMPTY = _dir_abritamr / 'html-empty.io'
OUTPUT_ABRITAMR_SUMMARY = _dir_abritamr / 'summary_out.tsv'
OUTPUT_ABRITAMR_SUMMARY_JSON = _dir_abritamr / 'summary_out.json'


def get_command_informs(config: Dict[str, Any]) -> List[Path]:
    """
    Returns the paths to the AbriTAMR informs.
    :param config: config
    :return: Path(s) to the AbriTAMR informs
    """
    paths = []

    # AbriTAMR is disabled -> return empty list
    if 'abritamr' not in config['analyses']:
        return []
    else:
        paths.append(OUTPUT_ABRITAMR_RUN_INFORMS)
        paths.append(OUTPUT_REPORT_ABRITAMR_INFORMS)

    return [Path(config['working_dir']) / p for p in paths]
