import numpy as np
from pathlib import Path
from typing import Any, List, Union

SNAKEFILE_GENOTYPHI = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_genotyphi = Path('genotyphi')
OUTPUT_GENOTYPHI = _dir_genotyphi / 'genotyphi_output.io'
OUTPUT_GENOTYPHI_INFORMS = _dir_genotyphi / 'informs.io'
OUTPUT_GENOTYPHI_REPORT = _dir_genotyphi / 'html.io'
OUTPUT_GENOTYPHI_REPORT_EMPTY = _dir_genotyphi / 'html-empty.io'
OUTPUT_GENOTYPHI_SUMMARY = _dir_genotyphi / 'summary_out.tsv'
OUTPUT_GENOTYPHI_SUMMARY_JSON = _dir_genotyphi / 'summary_out.json'


def numpy_to_python(np_type_value: Any) -> Union[int, float, str, List[Any]]:
    """
    Function to convert numpy types to python types.
    :param np_type_value: numpy type value
    :return: python type value
    """
    if isinstance(np_type_value, np.integer):
        return int(np_type_value)
    if isinstance(np_type_value, np.floating):
        return float(np_type_value)
    if isinstance(np_type_value, np.ndarray):
        return np_type_value.tolist()
    if isinstance(np_type_value, np.str_) or isinstance(np_type_value, str):
        return str(np_type_value)
