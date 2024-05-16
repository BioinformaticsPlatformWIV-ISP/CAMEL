from pathlib import Path
from typing import Dict, Any, List

from camel.app.loggers import logger

SNAKEFILE_VARIANT_CALLING = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_variant_calling = Path('variant_calling')

OUTPUT_VARIANT_CALLING_REPORT = _dir_variant_calling / 'report' / 'html.io'
OUTPUT_VARIANT_CALLING_REPORT_EMPTY = _dir_variant_calling / 'report' / 'html-empty.io'
OUTPUT_VARIANT_CALLING_SUMMARY = _dir_variant_calling / 'summary' / 'summary_out.tsv'
OUTPUT_VARIANT_CALLING_BAM = _dir_variant_calling / 'read_mapping' / 'bam.io'
OUTPUT_VARIANT_CALLING_DEPTH_INFORMS = _dir_variant_calling / 'depth' / 'informs.io'
OUTPUT_VARIANT_CALLING_DEPTH_TSV = _dir_variant_calling / 'depth' / 'tsv.io'
OUTPUT_VARIANT_CALLING_UNFILTERED_VCF = _dir_variant_calling / 'unzip_vcf' / 'vcf.io'
OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ = _dir_variant_calling / 'norm' / 'vcf_gz-indexed.io'
OUTPUT_VARIANT_CALLING_CONSENSUS = _dir_variant_calling / 'consensus' / 'fasta.io'
OUTPUT_VARIANT_CALLING_FILTERED_VCF = _dir_variant_calling / 'unzip_vcf_filtered' / 'vcf.io'
OUTPUT_VARIANT_CALLING_FILTERED_VCF_GZ = _dir_variant_calling / 'filter_zscore' / 'vcf_gz-indexed.io'
OUTPUT_VARIANT_CALLING_MAPPING_INFORMS = _dir_variant_calling / 'read_mapping' / 'informs.io'
OUTPUT_VARIANT_CALLING_INFORMS_ALL = _dir_variant_calling / 'informs_all.io'


def get_mapping_fq_input(config: Dict[str, Any]) -> Path:
    """
    Returns the fq input file path needed for mapping.
    :param config: Snakemake configuration
    :return: Path to the mapping fq input file
    """
    if config['input_type'] in ('illumina', 'ont', 'hybrid'): # Ont and hybrid were added because otherwise some tests of the mockpipeline fail
        return Path(config['working_dir']) / 'fq_dict.io'
    if config['input_type'] == 'fasta':
        return Path(config['working_dir']) / 'variant_calling' / 'art' / 'fastq.io'
    if config['input_type'] == 'fasta_with_vcf':  # was added because otherwise some tests of the mockpipeline fail
        logger.warning(f"Variant calling is not supported for input type '{config['input_type']}'")
        return Path('')


def get_bam(config: Dict[str, Any]) -> Path:
    """
    Returns the BAM output IO object path.
    :param config: Snakemake configuration
    :return: Path to the BAM file
    """
    if config['input_type'] in ('illumina', 'fasta', 'ont', 'hybrid'): # Ont and hybrid were added because otherwise some tests of the mockpipeline fail
        return Path(config['working_dir']) / OUTPUT_VARIANT_CALLING_BAM
    if config['input_type'] == 'fasta_with_vcf':
        return Path(config['working_dir']) / 'variant_calling' / 'dummy_bam' / 'bam.io'


def get_vcf(config: Dict[str, Any]) -> Path:
    """
    Returns the VCF output IO object path (before filtering).
    :param config: Snakemake configuration
    :return: Path to the unfiltered VCF file
    """
    if config['input_type'] in ('fasta', 'illumina', 'ont', 'hybrid'): # Ont and hybrid were added because otherwise some tests of the mockpipeline fail
        return OUTPUT_VARIANT_CALLING_UNFILTERED_VCF
    if config['input_type'] == 'fasta_with_vcf':
        return Path(config['working_dir']) / 'input' / 'vcf.io'


def get_vcf_gz(config: Dict[str, Any]) -> Path:
    """
    Returns the VCF GZ output IO object path (before filtering).
    :param config: Snakemake configuration
    :return: Path to the unfiltered gzipped VCF file
    """
    if config['input_type'] in ('fasta', 'illumina', 'ont', 'hybrid'): # Ont and hybrid were added because otherwise some tests of the mockpipeline fail
        return OUTPUT_VARIANT_CALLING_UNFILTERED_VCF_GZ
    if config['input_type'] == 'fasta_with_vcf':
        return Path(config['working_dir']) / 'variant_calling' / 'gzip' / 'vcf_gz.io'


def get_reports(config: Dict[str, Any]) -> Path:
    """
    Returns the path to the variant calling report.
    :param config: Snakemake configuration
    :return: Report path
    """
    input_type = config['input_type']

    if input_type in ('illumina', 'fasta'):
        return Path(config['working_dir']) / OUTPUT_VARIANT_CALLING_REPORT

    if input_type == 'fasta_with_vcf':
        return Path(config['working_dir']) / OUTPUT_VARIANT_CALLING_REPORT_EMPTY


def get_summaries(config: Dict[str, Any]) -> List[Path]:
    """
    Returns the paths to the variant calling summary file(s).
    :param config: Snakemake configuration
    :return: Summary path(s)
    """
    input_type = config['input_type']
    paths = []
    if input_type in ('illumina', 'fasta'):
        paths.append(OUTPUT_VARIANT_CALLING_SUMMARY)
    return [Path(config['working_dir']) / p for p in paths]


def get_command_informs(config: Dict[str, Any]) -> List[Path]:
    """
    Returns the paths to the variant calling informs.
    :param config: config
    :return: Path(s) to the variant calling informs
    """
    input_type = config['input_type']
    paths = []
    if input_type in ('illumina', 'fasta'):
        paths.append(OUTPUT_VARIANT_CALLING_INFORMS_ALL)
    return [Path(config['working_dir']) / p for p in paths]