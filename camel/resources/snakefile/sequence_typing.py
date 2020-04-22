from pathlib import Path
from typing import Dict, Any, Optional

SNAKEFILE_SEQUENCE_TYPING = f'{Path(__file__).parent / Path(__file__).stem}.smk'
SNAKEFILE_SEQUENCE_TYPING_BLAST = f'{Path(__file__).parent / Path(__file__).stem}_blast.smk'
SNAKEFILE_SEQUENCE_TYPING_SRST2 = f'{Path(__file__).parent / Path(__file__).stem}_srst2.smk'
SNAKEFILE_SEQUENCE_TYPING_KMA = f'{Path(__file__).parent / Path(__file__).stem}_kma.smk'

_dir_typing = Path('typing', '{scheme}')
OUTPUT_TYPING_REPORT = _dir_typing / 'html.io'
OUTPUT_TYPING_REPORT_EMPTY = _dir_typing / 'html-empty.io'
OUTPUT_TYPING_HITS = _dir_typing / '{locus_type}' / '{detection_method}' / 'all-hits.io'
OUTPUT_TYPING_SUMMARY = _dir_typing / 'summary_out.tsv'
OUTPUT_TYPING_INFORMS = _dir_typing / 'informs.io'


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
        return Path(config['working_dir']) / str(OUTPUT_TYPING_REPORT_EMPTY).format(scheme=scheme_key)
    return Path(config['working_dir']) / str(OUTPUT_TYPING_REPORT).format(scheme=scheme_key)
