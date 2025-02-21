from pathlib import Path
from typing import Any

from camel.resources.snakefile import read_simulation

SNAKEFILE_VARIANT_CALLING = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_variant_calling = Path('variant_calling')

OUTPUT_VARIANT_CALLING_REPORT = _dir_variant_calling / 'report' / 'html.io'
OUTPUT_VARIANT_CALLING_REPORT_EMPTY = _dir_variant_calling / 'report' / 'html-empty.io'
OUTPUT_VARIANT_CALLING_SUMMARY = _dir_variant_calling / 'summary' / 'summary_out.tsv'
OUTPUT_VARIANT_CALLING_DEPTH_INFORMS = _dir_variant_calling / 'depth' / 'informs.io'
OUTPUT_VARIANT_CALLING_DEPTH_TSV = _dir_variant_calling / 'depth' / 'tsv.io'
OUTPUT_VARIANT_CALLING_UNFILTERED_VCF = _dir_variant_calling / 'unzip_vcf' / 'vcf.io'
OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ = _dir_variant_calling / 'norm' / 'vcf_gz-indexed.io'
OUTPUT_VARIANT_CALLING_CONSENSUS = _dir_variant_calling / 'consensus' / 'fasta.io'
OUTPUT_VARIANT_CALLING_FILTERED_VCF = _dir_variant_calling / 'unzip_vcf_filtered' / 'vcf.io'
OUTPUT_VARIANT_CALLING_FILTERED_VCF_GZ = _dir_variant_calling / 'filter_zscore' / 'vcf_gz-indexed.io'
OUTPUT_VARIANT_CALLING_MAPPING_RATE_INFORMS = _dir_variant_calling / 'rate' / 'informs.io'
OUTPUT_VARIANT_CALLING_MAPPING_INFORMS = _dir_variant_calling / 'read_mapping' / 'informs.io'
OUTPUT_VARIANT_CALLING_INFORMS_ALL = _dir_variant_calling / 'informs_all.io'


def get_mapping_fq_input(config: dict[str, Any]) -> Path:
    """
    Returns the fq input file path needed for mapping.
    :param config: Snakemake configuration
    :return: Path to the mapping fq input file
    """
    # ONT and hybrid were added because otherwise some tests of the MockPipeline fail
    if config['input_type'] in ('illumina', 'ont', 'hybrid'):
        return Path(config['working_dir']) / 'fq_dict.io'
    if config['input_type'] in ('fasta', 'fasta_with_vcf'):
        return Path(config['working_dir']) / read_simulation.OUTPUT_SIMULATION_FASTQ


def get_bam(config: dict[str, Any]) -> Path:
    """
    Returns the BAM output IO object path.
    :param config: Snakemake configuration
    :return: Path to the BAM file
    """
    # ONT and hybrid were added because otherwise some tests of the MockPipeline fail
    if config['input_type'] in ('illumina', 'fasta', 'hybrid'):
        return Path(config['working_dir']) / 'variant_calling' / 'read_mapping' / 'illumina' / 'bam.io'  # OUTPUT_VARIANT_CALLING_BAM_ILLUMINA
    if config['input_type'] == 'ont':
        return Path(config['working_dir']) / 'variant_calling' / 'read_mapping' / 'ont' / 'bam.io'  # OUTPUT_VARIANT_CALLING_BAM_ONT
    if config['input_type'] == 'fasta_with_vcf':
        return Path(config['working_dir']) / 'variant_calling' / 'dummy_bam' / 'bam.io'


def get_vcf(config: dict[str, Any]) -> Path:
    """
    Returns the VCF output IO object path (before filtering).
    :param config: Snakemake configuration
    :return: Path to the unfiltered VCF file
    """
    if config['input_type'] in ('fasta', 'illumina', 'ont', 'hybrid'):  # Ont and hybrid were added because otherwise some tests of the mockpipeline fail
        return OUTPUT_VARIANT_CALLING_UNFILTERED_VCF
    if config['input_type'] == 'fasta_with_vcf':
        return Path(config['working_dir']) / 'input' / 'vcf.io'


def get_vcf_gz(config: dict[str, Any]) -> Path:
    """
    Returns the VCF GZ output IO object path (before filtering).
    :param config: Snakemake configuration
    :return: Path to the unfiltered gzipped VCF file
    """
    if config['input_type'] in ('fasta', 'illumina', 'ont', 'hybrid'):  # Ont and hybrid were added because otherwise some tests of the mockpipeline fail
        return OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ
    if config['input_type'] == 'fasta_with_vcf':
        return Path(config['working_dir']) / 'variant_calling' / 'gzip' / 'vcf_gz.io'


def get_reports(config: dict[str, Any]) -> Path:
    """
    Returns the path to the variant calling report.
    :param config: Snakemake configuration
    :return: Report path
    """
    input_type = config['input_type']

    if input_type in ('illumina', 'fasta', 'ont'):
        return Path(config['working_dir']) / OUTPUT_VARIANT_CALLING_REPORT

    if input_type == 'fasta_with_vcf':
        return Path(config['working_dir']) / OUTPUT_VARIANT_CALLING_REPORT_EMPTY


def get_mapping_informs(config: dict[str, Any]) -> Path:
    """
    Returns the paths to the variant calling mapping informs.
    :param config: config
    :return: Path(s) to the variant calling mapping informs
    """
    input_type = config['input_type']

    if input_type in ('illumina', 'fasta', 'fasta_with_vcf'):
        return Path(config['working_dir']) / 'variant_calling' / 'read_mapping' / 'illumina' / 'informs.io'
    elif input_type == 'ont':
        return Path(config['working_dir']) / 'variant_calling' / 'read_mapping' / 'ont' / 'informs.io'
    else:
        raise ValueError(f"No read mapping for input {input_type}")

def get_summaries(config: dict[str, Any]) -> list[Path]:
    """
    Returns the paths to the variant calling summary file(s).
    :param config: Snakemake configuration
    :return: Summary path(s)
    """
    input_type = config['input_type']
    paths = []
    if input_type in ('illumina', 'fasta', 'ont'):
        paths.append(OUTPUT_VARIANT_CALLING_SUMMARY)
    return [Path(config['working_dir']) / p for p in paths]


def get_command_informs(config: dict[str, Any]) -> list[Path]:
    """
    Returns the paths to the variant calling informs.
    :param config: config
    :return: Path(s) to the variant calling informs
    """
    input_type = config['input_type']
    paths = []
    if input_type in ('illumina', 'fasta', 'ont'):
        paths.append(OUTPUT_VARIANT_CALLING_INFORMS_ALL)
    return [Path(config['working_dir']) / p for p in paths]
