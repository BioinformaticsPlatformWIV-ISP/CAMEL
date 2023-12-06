from pathlib import Path
from typing import List

SNAKEFILE_CONFINDR = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_confindr = Path('confindr')
OUTPUT_CONFINDR_INFORMS = _dir_confindr / 'informs.io'
OUTPUT_CONFINDR_REPORT = _dir_confindr / 'report' / 'html.io'
OUTPUT_CONFINDR_REPORT_EMPTY = _dir_confindr / 'report' / 'html-empty.io'
OUTPUT_CONFINDR_SUMMARY = _dir_confindr / 'summary' / 'summary_confindr.tsv'


def get_report(config) -> Path:
    """
    Returns the path to the ConFindr report io file.
    :param config: Snakemake configuration
    :return: Path to report
    """
    if ('confindr' not in config['analyses']) or (config['input_type'] not in ['illumina', 'hybrid']):
        return Path(config['working_dir'], OUTPUT_CONFINDR_REPORT_EMPTY)
    return Path(config['working_dir'], OUTPUT_CONFINDR_REPORT)


def get_command_informs(config) -> List[Path]:
    """
    Returns the path to the ConFindr informs io file.
    :param config: Snakemake configuration
    :return: Path to informs IO object
    """
    if ('confindr' not in config['analyses']) or (config['input_type'] not in ['illumina', 'hybrid']):
        return []
    return [Path(config['working_dir'], OUTPUT_CONFINDR_INFORMS)]


def get_summary(config) -> List[Path]:
    """
    Returns the path to the ConFindr summay file.
    :param config: Snakemake configuration
    :return: Path to summary TSV
    """
    if ('confindr' not in config['analyses']) or (config['input_type'] not in ['illumina', 'hybrid']):
        return []
    return [Path(config['working_dir'], OUTPUT_CONFINDR_SUMMARY)]
