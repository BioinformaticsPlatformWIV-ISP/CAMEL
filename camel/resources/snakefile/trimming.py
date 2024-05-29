from pathlib import Path
from typing import Dict, Any, List

from camel.resources.snakefile import trimming_illumina, trimming_ont


def get_reports(config: Dict[str, Any]) -> List[Path]:
    """
    Returns the path to the read trimming report.
    :param config: Snakemake configuration
    :return: List of paths to the trimming reports
    """
    paths = []
    input_type = config['input_type']

    # FASTA input
    if input_type in ('fasta', 'fasta_with_vcf'):
        return []

    # FASTQ input
    if input_type in ('illumina', 'hybrid'):
        paths.append(trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_REPORT)
    if input_type in ('ont', 'hybrid'):
        paths.append(trimming_ont.OUTPUT_TRIMMING_ONT_REPORT)

    # Check if at least one path was added
    if len(paths) == 0:
        raise ValueError(f'Invalid input type: {input_type}')

    return [Path(config['working_dir']) / p for p in paths]


def get_summaries(config: Dict[str, Any]) -> List[Path]:
    """
    Returns the paths to the read trimming summary file(s).
    :param config: Snakemake configuration
    :return: Path to read trimming summary file
    """
    input_type = config['input_type']

    # FASTA input
    if input_type in ('fasta', 'fasta_with_vcf'):
        return []

    # FASTQ input
    paths = []
    if input_type in ('illumina', 'hybrid'):
        paths.append(trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_SUMMARY)
    if input_type in ('ont', 'hybrid'):
        paths.append(trimming_ont.OUTPUT_TRIMMING_ONT_SUMMARY)
    return [Path(config['working_dir']) / p for p in paths]


def get_command_informs(config: Dict[str, Any]) -> List[Path]:
    """
    Returns a list of informs IO files for the read trimming steps.
    :param config: Snakemake configuration
    :return: List of informs IO files
    """
    paths = []
    input_type = config['input_type']
    if input_type in ('illumina', 'hybrid'):
        paths.append(trimming_illumina.OUTPUT_TRIMMING_ILLUMINA_INFORMS)
    if input_type in ('ont', 'hybrid'):
        paths.append(trimming_ont.OUTPUT_TRIMMING_ONT_INFORMS)
    return [Path(config['working_dir']) / p for p in paths]
