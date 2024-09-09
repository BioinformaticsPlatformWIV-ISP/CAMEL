from pathlib import Path
from typing import Dict, Any, List

SNAKEFILE_SCRUBBING = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_scrubbing = Path('human_read_scrubbing', '{input_format}')
INPUT_SCRUBBING_FASTQ = _dir_scrubbing / 'input' / 'fastq.io'
INPUT_SCRUBBING_FASTA = _dir_scrubbing / 'input' / 'fasta.io'
OUTPUT_SCRUBBING_INFORMS = _dir_scrubbing / 'scrubbing' / 'informs.io'
OUTPUT_SCRUBBING_REPORT = _dir_scrubbing / 'output' / 'html.io'
OUTPUT_SCRUBBING_REPORT_EMPTY = _dir_scrubbing / 'output' / 'html-empty.io'
OUTPUT_SCRUBBING_SUMMARY = _dir_scrubbing / 'output' / 'summary_out.tsv'
OUTPUT_SCRUBBING_SUMMARY_JSON = _dir_scrubbing / 'output' / 'summary_out.json'
OUTPUT_SCRUBBING_FASTQ = _dir_scrubbing / 'output' / 'fastq.io'
OUTPUT_SCRUBBING_FASTA = _dir_scrubbing / 'output' / 'fasta.io'


def get_reports(config: Dict[str, Any]) -> List[Path]:
    """
    Returns the paths to the human read scrubbing reports.
    :param config: Snakemake configuration
    :return: Report path(s)
    """
    input_type = config['input_type']
    paths = []

    # FASTA input
    if input_type in ('fasta', 'fasta_with_vcf'):
        if 'human_read_scrubbing' in config['analyses']:
            paths.append(str(OUTPUT_SCRUBBING_REPORT).format(input_format='fasta'))
        else:
            paths.append(str(OUTPUT_SCRUBBING_REPORT_EMPTY).format(input_format='fasta'))

    # PE reads
    if input_type in ('illumina', 'hybrid'):
        if 'human_read_scrubbing' in config['analyses']:
            paths.append(str(OUTPUT_SCRUBBING_REPORT).format(input_format='fastq_pe'))
        else:
            paths.append(str(OUTPUT_SCRUBBING_REPORT_EMPTY).format(input_format='fastq_pe'))

    # SE reads
    if input_type in ('ont', 'hybrid'):
        if 'human_read_scrubbing' in config['analyses']:
            paths.append(str(OUTPUT_SCRUBBING_REPORT).format(input_format='fastq_se'))
        else:
            paths.append(str(OUTPUT_SCRUBBING_REPORT_EMPTY).format(input_format='fastq_se'))

    return [Path(config['working_dir']) / p for p in paths]


def get_command_informs(config: Dict[str, Any]) -> List[Path]:
    """
    Returns a list of informs IO files for the human read scrubbing steps.
    :param config: Snakemake configuration
    :return: List of informs IO files
    """
    input_type = config['input_type']
    paths = []

    if 'human_read_scrubbing' not in config['analyses']:
        return []

    # FASTA input
    if (input_type in ('fasta', 'fasta_with_vcf')) and ('human_read_scrubbing' in config['analyses']):
        paths.append(str(OUTPUT_SCRUBBING_INFORMS).format(input_format='fasta'))

    # PE reads
    if (input_type in ('illumina', 'hybrid')) and ('human_read_scrubbing' in config['analyses']):
        paths.append(str(OUTPUT_SCRUBBING_INFORMS).format(input_format='fastq_pe'))

    # SE reads
    if (input_type in ('ont', 'hybrid')) and ('human_read_scrubbing' in config['analyses']):
        paths.append(str(OUTPUT_SCRUBBING_INFORMS).format(input_format='fastq_se'))
    return [Path(config['working_dir']) / p for p in paths]


def get_summaries(config: Dict[str, Any]) -> List[Path]:
    """
    Returns the paths to the human read scrubbing summary file(s).
    :param config: Snakemake configuration
    :return: Summary file path(s)
    """
    input_type = config['input_type']
    paths = []

    if 'human_read_scrubbing' not in config['analyses']:
        return []

    # FASTA input
    if (input_type in ('fasta', 'fasta_with_vcf')) and ('human_read_scrubbing' in config['analyses']):
        paths.append(str(OUTPUT_SCRUBBING_SUMMARY).format(input_format='fasta'))

    # PE reads
    if (input_type in ('illumina', 'hybrid')) and ('human_read_scrubbing' in config['analyses']):
        paths.append(str(OUTPUT_SCRUBBING_SUMMARY).format(input_format='fastq_pe'))

    # SE reads
    if (input_type in ('ont', 'hybrid')) and ('human_read_scrubbing' in config['analyses']):
        paths.append(str(OUTPUT_SCRUBBING_SUMMARY).format(input_format='fastq_se'))

    return [Path(config['working_dir']) / p for p in paths]


def get_output_io(config: Dict[str, Any]) -> Path:
    """
    Returns the paths to the human read scrubbing output io file(s).
    :param config: Snakemake configuration
    :return: Summary file path(s)
    """
    input_type = config['input_type']

    # FASTA input
    if (input_type in ('fasta', 'fasta_with_vcf')) and ('human_read_scrubbing' in config['analyses']):
        return Path(config['working_dir']) / str(OUTPUT_SCRUBBING_FASTA).format(input_format='fasta')

    # PE reads
    if (input_type in ('illumina', 'hybrid')) and ('human_read_scrubbing' in config['analyses']):
        return Path(config['working_dir']) / str(OUTPUT_SCRUBBING_FASTQ).format(input_format='fastq_pe')

    # SE reads
    if (input_type in ('ont', 'hybrid')) and ('human_read_scrubbing' in config['analyses']):
        return Path(config['working_dir']) / str(OUTPUT_SCRUBBING_FASTQ).format(input_format='fastq_se')
