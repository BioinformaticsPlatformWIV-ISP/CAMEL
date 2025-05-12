import re
from pathlib import Path
from typing import List

SNAKEFILE_REF_SELECTION = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_ref_selection = Path('ref_selection')

OUTPUT_REF_SELECTION_REPORT = _dir_ref_selection / 'report' / 'html.io'
OUTPUT_REF_SELECTION_REPORT_EMPTY = _dir_ref_selection / 'report' / 'html-empty.io'
OUTPUT_REF_SELECTION_FASTA = _dir_ref_selection / 'create_fasta' / 'fasta.io'
OUTPUT_REF_SELECTION_SUMMARY = _dir_ref_selection / 'summary_ref_selection.tsv'


def get_segments(path_db: Path) -> List[str]:
    """
    Returns the names of the segments in the current database.
    :param path_db: Path to the database
    """
    segments = []
    for file_ in (path_db / 'mash').iterdir():
        if not file_.name.endswith('.msh'):
            continue
        segments.append(re.search(r'(.*)\.msh', file_.name).group(1))
    return sorted(segments)
