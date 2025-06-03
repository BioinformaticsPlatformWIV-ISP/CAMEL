from pathlib import Path
from typing import Dict, Any, List

SNAKEFILE_ABRITAMR = Path(__file__).parent / f'{Path(__file__).stem}.smk'
_dir_abritamr = Path('abritamr')
OUTPUT_ABRITAMR_MATCHES = _dir_abritamr / 'abritamr_output_matches.io'
OUTPUT_ABRITAMR_PARTIALS = _dir_abritamr / 'abritamr_output_partials.io'
OUTPUT_ABRITAMR_COMBINED = _dir_abritamr / 'abritamr_output_combined.io'
OUTPUT_ABRITAMR_QC = _dir_abritamr / 'qc_file.io'
OUTPUT_ABRITAMR_RUN_INFORMS = _dir_abritamr / 'informs_run.io'
OUTPUT_ABRITAMR_REPORT_REPORT = _dir_abritamr / 'abritamr_output_report.io'
OUTPUT_ABRITAMR_REPORT_REPORT_INFORMS = _dir_abritamr / 'informs_report.io'
OUTPUT_ABRITAMR_REPORT = _dir_abritamr / 'html.io'
OUTPUT_ABRITAMR_REPORT_EMPTY = _dir_abritamr / 'html-empty.io'
OUTPUT_ABRITAMR_SUMMARY = _dir_abritamr / 'summary_out.tsv'


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
        paths.append(OUTPUT_ABRITAMR_REPORT_REPORT_INFORMS)

    return [Path(config['working_dir']) / p for p in paths]
