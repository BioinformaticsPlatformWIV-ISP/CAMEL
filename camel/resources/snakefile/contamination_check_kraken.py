from pathlib import Path
from typing import Dict

SNAKEFILE_CONTAMINATION_CHECK_KRAKEN = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_out = Path('contamination_check', '{read_type}')
OUTPUT_CONTAMINATION_CHECK_REPORT = _dir_out / 'report' / 'html.io'
OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY = _dir_out / 'report' / 'html-empty.io'
OUTPUT_CONTAMINATION_CHECK_INFORMS = _dir_out / 'kraken2' / 'informs-contamination.io'
OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS = _dir_out / 'kraken2' / 'informs.io'
OUTPUT_CONTAMINATION_SUMMARY = _dir_out / 'summary' / 'summary_out.tsv'


def get_contamination_check_report(config, read_type: str = None) -> Path:
    """
    Returns the Kraken report.
    :param config: config
    :param read_type: read_type
    :return: Path to the Kraken informs
    """
    if read_type is None:
        read_type = config['read_type']
    if read_type not in ['illumina', 'nanopore']:
        raise ValueError(f'Invalid read type: {read_type}')
    if 'kraken' not in config['analyses']:
        return Path(config['working_dir']) / str(OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY).format(read_type=read_type)
    else:
        return Path(config['working_dir']) / str(OUTPUT_CONTAMINATION_CHECK_REPORT).format(read_type=read_type)


def get_contamination_check_summary(config: Dict, read_type: str = None) -> Path:
    """
    Returns the Kraken summary.
    :param config: config
    :param read_type: read_type
    :return: Path to the Kraken informs
    """
    if read_type is None:
        read_type = config['read_type']
    if read_type not in ['illumina', 'nanopore']:
        raise ValueError(f'Invalid read type: {read_type}')
    return Path(config['working_dir']) / str(OUTPUT_CONTAMINATION_SUMMARY).format(read_type=read_type)


def get_contamination_check_kraken_informs(config: Dict, read_type: str = None) -> Path:
    """
    Returns the Kraken informs.
    :param config: config
    :param read_type: read_type
    :return: Path to the Kraken informs
    """
    if read_type is None:
        read_type = config['read_type']
    if read_type not in ['illumina', 'nanopore']:
        raise ValueError(f'Invalid read type: {read_type}')
    return Path(config['working_dir']) / str(OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS).format(read_type=read_type)
