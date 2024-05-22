from pathlib import Path
from typing import Dict, Any

SNAKEFILE_SHIGATYPER = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_shigatyper = Path('shigatyper')
OUTPUT_SHIGATYPER_REPORT = _dir_shigatyper/ 'report' / 'html.io'
OUTPUT_SHIGATYPER_REPORT_EMPTY = _dir_shigatyper / 'report' / 'html-empty.io'
OUTPUT_SHIGATYPER_SUMMARY = _dir_shigatyper / 'summary_shigatyper.tsv'
OUTPUT_SHIGATYPER_INFORMS = _dir_shigatyper/ 'informs.io'


def get_input(config: Dict[str, Any]) -> Path:
    """
    Returns the input for the ShigaTyper tool.
    :param config: Snakemake configuration
    :return: Path to the input file
    """
    if config['input_type'] == 'illumina':
        return Path(config['working_dir']) / 'fq_dict.io'
    if config['input_type'] == 'fasta':
        return Path(config['working_dir']) / 'variant_calling' / 'art' / 'fastq.io'
    