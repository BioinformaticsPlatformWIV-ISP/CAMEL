from pathlib import Path
from typing import Dict, Any, List

SNAKEFILE_DOWNSAMPLING = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_downsampling = Path('downsampling', '{read_key}')
INPUT_DOWNSAMPLING_FASTQ = _dir_downsampling / 'input' / 'fastq.io'
OUTPUT_DOWNSAMPLING_REPORT = _dir_downsampling / 'html.io'
OUTPUT_DOWNSAMPLING_SUMMARY = _dir_downsampling / 'summary_downsampling.tsv'
OUTPUT_DOWNSAMPLING_FASTQ = _dir_downsampling / 'seqtk' / 'fastq.io'
OUTPUT_DOWNSAMPLING_INFORMS = _dir_downsampling / 'seqtk' / 'informs.io'


def get_reports(config: Dict[str, Any], input_type: str = None) -> List[Path]:
    """
    Returns the paths to the downsampling reports.
    :param input_type: input type
    :param config: Snakemake configuration
    :return: Report path(s)
    """
    if input_type is None:
        input_type = config['input_type']

    # FASTA input
    if input_type in ('fasta', 'fasta_with_vcf'):
        return []

    # FASTQ input
    paths = []
    if input_type in ('illumina', 'hybrid'):
        paths.append(str(OUTPUT_DOWNSAMPLING_REPORT).format(read_key='fastq_pe'))
    if input_type in ('ont', 'hybrid', 'iontorrent'):
        paths.append(str(OUTPUT_DOWNSAMPLING_REPORT).format(read_key='fastq_se'))

    # Check if at least one path was added
    if len(paths) == 0:
        raise ValueError(f'No downsampling report(s) for input type: {input_type}')

    return [Path(config['working_dir']) / p for p in paths]


def get_command_informs(config: Dict[str, Any], input_type: str = None) -> List[Path]:
    """
    Returns a list of informs IO files for the downsampling steps.
    :param config: Snakemake configuration
    :param input_type: Input type
    :return: List of informs IO files
    """
    if input_type is None:
        input_type = config['input_type']

    # FASTA input
    if input_type in ('fasta', 'fasta_with_vcf'):
        return []

    # FASTQ input
    paths = []
    if input_type in ('hybrid', 'illumina'):
        paths.append(str(OUTPUT_DOWNSAMPLING_INFORMS).format(read_key='fastq_pe'))
    if input_type in ('iontorrent', 'hybrid', 'ont'):
        paths.append(str(OUTPUT_DOWNSAMPLING_INFORMS).format(read_key='fastq_se'))

    # Check if at least one path was added
    if len(paths) == 0:
        raise ValueError(f'No downsampling informs for input type: {input_type}')

    return [Path(config['working_dir']) / p for p in paths]


def get_summaries(config: Dict[str, Any], input_type: str = None) -> List[Path]:
    """
    Returns the paths to the downsampling summary file(s).
    :param input_type: input type
    :param config: Snakemake configuration
    :return: Summary file path(s)
    """
    if input_type is None:
        input_type = config['input_type']

    # FASTA input
    if input_type in ('fasta', 'fasta_with_vcf'):
        return []

    # FASTQ input
    paths = []
    if input_type in ('hybrid', 'illumina'):
        paths.append(str(OUTPUT_DOWNSAMPLING_SUMMARY).format(read_key='fastq_pe'))
    if input_type in ('iontorrent', 'hybrid', 'ont'):
        paths.append(str(OUTPUT_DOWNSAMPLING_SUMMARY).format(read_key='fastq_se'))

    # Check if at least one path was added
    if len(paths) == 0:
        raise ValueError(f'No downsampling summary output for input type: {input_type}')

    return [Path(config['working_dir']) / p for p in paths]
