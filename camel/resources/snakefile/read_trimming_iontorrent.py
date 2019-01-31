from typing import Dict, Any

import os

from camel.resources.snakefile.read_trimming import OUTPUT_READ_TRIMMING_REPORT, FOLDER_TRIMMING

FOLDER_TRIMMING_IT = 'read_trimming_it'
OUTPUT_TRIMMING_IT_REPORT = os.path.join('read_trimming_it', 'report', 'html.io')
OUTPUT_TRIMMING_IT_SUMMARY = os.path.join('read_trimming_it', 'summary', 'summary_out.tsv')
OUTPUT_TRIMMING_IT_READS = os.path.join('read_trimming_it', 'trim_qual', 'fastq.io')


def get_read_trimming_report(config: Dict[str, Any]) -> str:
    """
    Returns the path to the read trimming report.
    :param config: Snakemake configuration
    :return: Path to read trimming report
    """
    if ('read_type' not in config) or (config['read_type'] == 'illumina'):
        relative_path = OUTPUT_READ_TRIMMING_REPORT
    elif config['read_type'] == 'iontorrent':
        relative_path = OUTPUT_TRIMMING_IT_REPORT
    else:
        raise ValueError(f"Invalid read type: '{config['read_type']}'")
    return os.path.join(config['working_dir'], relative_path)


def get_trimming_folder(config: Dict[str, Any]) -> str:
    """
    Returns the read trimming folder.
    :param config: Snakemake configuration
    :return: Path to read trimming folder
    """
    if ('read_type' not in config) or (config['read_type'] == 'illumina'):
        relative_path = FOLDER_TRIMMING
    elif config['read_type'] == 'iontorrent':
        relative_path = FOLDER_TRIMMING_IT
    else:
        raise ValueError(f"Invalid read type: '{config['read_type']}'")
    return os.path.join(config['working_dir'], relative_path)
