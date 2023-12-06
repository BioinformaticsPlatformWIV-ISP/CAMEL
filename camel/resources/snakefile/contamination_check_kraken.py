from pathlib import Path
from typing import Dict, List, Any

SNAKEFILE_CONTAMINATION_CHECK_KRAKEN = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_out = Path('contamination_check', '{read_key}')
OUTPUT_CONTAMINATION_CHECK_REPORT = _dir_out / 'report' / 'html.io'
OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY = _dir_out / 'report' / 'html-empty.io'
OUTPUT_CONTAMINATION_CHECK_INFORMS = _dir_out / 'kraken2' / 'informs-contamination.io'
OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS = _dir_out / 'kraken2' / 'informs.io'
OUTPUT_CONTAMINATION_SUMMARY = _dir_out / 'summary' / 'summary_out.tsv'


def get_reports(config: Dict[str, Any]) -> List[Path]:
    """
    Returns the paths to the contamination check reports.
    :param config: Snakemake configuration
    :return: Report path(s)
    """
    input_type = config['input_type']
    paths = []

    # PE reads
    if input_type in ('illumina', 'hybrid'):
        if 'kraken2' in config['analyses']:
            paths.append(str(OUTPUT_CONTAMINATION_CHECK_REPORT).format(read_key='fastq_pe'))
        else:
            paths.append(str(OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY).format(read_key='fastq_pe'))

    # SE reads
    if input_type in ('ont', 'hybrid'):
        if 'kraken2' in config['analyses']:
            paths.append(str(OUTPUT_CONTAMINATION_CHECK_REPORT).format(read_key='fastq_se'))
        else:
            paths.append(str(OUTPUT_CONTAMINATION_CHECK_REPORT_EMPTY).format(read_key='fastq_se'))
    return [Path(config['working_dir']) / p for p in paths]


def get_summaries(config: Dict[str, Any]) -> List[Path]:
    """
    Returns the paths to the contamination check summary file(s).
    :param config: Snakemake configuration
    :return: Report path(s)
    """
    input_type = config['input_type']
    paths = []

    # PE reads
    if (input_type in ('illumina', 'hybrid')) and ('kraken2' in config['analyses']):
        paths.append(str(OUTPUT_CONTAMINATION_SUMMARY).format(read_key='fastq_pe'))

    # SE reads
    if (input_type in ('ont', 'hybrid')) and ('kraken2' in config['analyses']):
        paths.append(str(OUTPUT_CONTAMINATION_SUMMARY).format(read_key='fastq_se'))
    return [Path(config['working_dir']) / p for p in paths]


def get_command_informs(config: Dict[str, Any]) -> List[Path]:
    """
    Returns the paths to the Kraken informs.
    :param config: config
    :return: Path(s) to the Kraken informs
    """
    input_type = config['input_type']
    paths = []

    # Kraken 2 is disabled -> return empty list
    if 'kraken2' not in config['analyses']:
        return []

    # PE reads
    if input_type in ('illumina', 'hybrid'):
        paths.append(str(OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS).format(read_key='fastq_pe'))

    # SE reads
    if input_type in ('ont', 'hybrid'):
        paths.append(str(OUTPUT_CONTAMINATION_CHECK_KRAKEN_INFORMS).format(read_key='fastq_se'))
    return [Path(config['working_dir']) / p for p in paths]
