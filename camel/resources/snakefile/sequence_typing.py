from typing import Dict, Any, Optional

import os

OUTPUT_TYPING_REPORT = os.path.join('typing', '{scheme}', 'html.io')
OUTPUT_TYPING_HITS = os.path.join('typing', '{scheme}', '{locus_type}', '{detection_method}', 'all-hits.io')
OUTPUT_TYPING_REPORT_EMPTY = os.path.join('typing', '{scheme}', 'html-empty.io')
OUTPUT_TYPING_SUMMARY = os.path.join('typing', '{scheme}', 'summary_out.tsv')


def get_sequence_typing_report(scheme_key: str, config: Dict[str, Any], analysis_name: Optional[str] = None) -> str:
    """
    Returns the report input for the given database key.
    :param scheme_key: Database key
    :param config: Pipeline config
    :param analysis_name: Analysis name that is checked
    :return: Report input path
    """
    search_key = analysis_name if analysis_name is not None else scheme_key
    if search_key not in config['analyses']:
        return os.path.join(config['working_dir'], OUTPUT_TYPING_REPORT_EMPTY.format(scheme=scheme_key))
    return os.path.join(config['working_dir'], OUTPUT_TYPING_REPORT.format(scheme=scheme_key))
