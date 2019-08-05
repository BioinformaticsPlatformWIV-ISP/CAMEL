from typing import Dict, Any, List

import os

from camel.resources.snakefile.read_trimming import OUTPUT_READ_TRIMMING_REPORT, FOLDER_TRIMMING, \
    OUTPUT_READ_TRIMMING_SUMMARY, OUTPUT_READ_TRIMMING_INFORMS

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


def get_read_trimming_summary(config: Dict[str, Any]) -> str:
    """
    Returns the path to the read trimming summary file.
    :param config: Snakemake configuration
    :return: Path to read trimming summary file
    """
    if ('read_type' not in config) or (config['read_type'] == 'illumina'):
        relative_path = OUTPUT_READ_TRIMMING_SUMMARY
    elif config['read_type'] == 'iontorrent':
        relative_path = OUTPUT_TRIMMING_IT_SUMMARY
    else:
        raise ValueError(f"Invalid read type: '{config['read_type']}'")
    return os.path.join(config['working_dir'], relative_path)


def get_read_trimming_commands(config: Dict[str, Any]) -> List[str]:
    """
    Returns a list of informs IO files for the read trimming steps.
    :return: List of informs IO files
    """
    if ('read_type' not in config) or (config['read_type'] == 'illumina'):
        return [os.path.join(config['working_dir'], OUTPUT_READ_TRIMMING_INFORMS)]
    elif config['read_type'] == 'iontorrent':
        return [
            os.path.join(config['working_dir'], 'read_trimming_it', 'trim_length', 'informs.io'),
            os.path.join(config['working_dir'], 'read_trimming_it', 'trim_qual', 'informs.io')
        ]
    else:
        raise ValueError(f"Invalid read type: '{config['read_type']}'")


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
