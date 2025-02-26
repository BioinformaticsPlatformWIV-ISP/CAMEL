from pathlib import Path
from typing import Dict, Any

from camel.resources.snakefile import assembly

SNAKEFILE_MYKROBE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_mykrobe = Path('mykrobe')
OUTPUT_MYKROBE_REPORT = _dir_mykrobe / 'report' / 'html.io'
OUTPUT_MYKROBE_REPORT_EMPTY = _dir_mykrobe / 'report' / 'html-empty.io'
OUTPUT_MYKROBE_SUMMARY = _dir_mykrobe / 'summary_mykrobe.tsv'
OUTPUT_MYKROBE_INFORMS = _dir_mykrobe / 'informs.io'


def get_input(config: Dict[str, Any]) -> Path:
    """
    Returns the input for the Mykrobe tool.
    :param config: Snakemake configuration
    :return: Path to the input file
    """
    if config['input_type'] in ('illumina', 'ont'):
        return Path(config['working_dir']) / 'fq_dict.io'
    if config['input_type'] == 'fasta':
        return Path(config['working_dir']) / assembly.OUTPUT_ASSEMBLY_FASTA
