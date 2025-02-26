from pathlib import Path
from typing import Dict, Any, List

from camel.app.loggers import logger
from camel.resources.snakefile import assembly_spades, assembly_flye, polish_assembly_short, polish_assembly_long, \
    human_read_scrubbing
from camel.resources.snakefile.variant_calling import OUTPUT_VARIANT_CALLING_MAPPING_INFORMS, \
    OUTPUT_VARIANT_CALLING_DEPTH_INFORMS

SNAKEFILE_ASSEMBLY = f'{Path(__file__).parent / Path(__file__).stem}.smk'
OUTPUT_ASSEMBLY_FASTA = Path('assembly', 'filtering', 'fasta.io')
OUTPUT_ASSEMBLY_FILTERING_INFORMS = Path('assembly', 'filtering', 'informs.io')
OUTPUT_ASSEMBLY_REPORT = Path('assembly', 'report', 'html.io')
OUTPUT_ASSEMBLY_SUMMARY = Path('assembly', 'summary', 'summary_out.tsv')
OUTPUT_ASSEMBLY_MAPPING_INFORMS = Path('assembly', '{mapper}', 'informs.io')
OUTPUT_ASSEMBLY_DEPTH_INFORMS = Path('assembly', 'samtools_depth', '{mapper}', 'informs.io')
OUTPUT_ASSEMBLY_MAPPING_RATE_INFORMS = Path('assembly', 'samtools_flagstat', '{mapper}', 'informs.io')


def get_fasta_raw(config: Dict[str, Any]) -> Path:
    """
    Returns the assembly FASTA output IO object path (before filtering).
    """
    if (config['input_type'] in ('fasta', 'fasta_with_vcf')) and ('human_read_scrubbing' not in config['analyses']):
        return Path(str(human_read_scrubbing.INPUT_SCRUBBING_FASTA).format(input_format='fasta'))
    if (config['input_type'] in ('fasta', 'fasta_with_vcf')) and ('human_read_scrubbing' in config['analyses']):
        return Path(str(human_read_scrubbing.OUTPUT_SCRUBBING_FASTA).format(input_format='fasta'))
    if config['input_type'] == 'illumina':
        return assembly_spades.OUTPUT_ASSEMBLY_FASTA
    if config['input_type'] == 'ont':
        return assembly_flye.OUTPUT_ASSEMBLY_FASTA
    if config['input_type'] == 'hybrid':
        return Path(str(polish_assembly_short.OUTPUT_POLISHING_FASTA).format(assembly_type='flye'))
    raise ValueError(f"Invalid input type: {config['input_type']}")


def get_command_informs(config: Dict[str, Any]) -> List[Path]:
    """
    Returns the assembly informs output IO object paths.
    :return: Assembly informs path
    """
    if config['input_type'] in ('fasta', 'fasta_with_vcf'):
        return []
    if config['input_type'] == 'illumina':
        return [Path(config['working_dir'], assembly_spades.OUTPUT_ASSEMBLY_INFORMS)]
    if config['input_type'] == 'ont':
        return [Path(config['working_dir'], assembly_flye.OUTPUT_ASSEMBLY_INFORMS)]
    if config['input_type'] == 'hybrid':
        return [
            Path(config['working_dir'], assembly_flye.OUTPUT_ASSEMBLY_INFORMS),
            Path(config['working_dir'], str(polish_assembly_long.OUTPUT_POLISH_MEDAKA_INFORMS).format(assembly_type='flye')),
            Path(config['working_dir'], str(polish_assembly_short.OUTPUT_POLYPOLISH_INFORMS).format(assembly_type='flye')),
            Path(config['working_dir'], str(polish_assembly_short.OUTPUT_POLCA_INFORMS).format(assembly_type='flye'))
        ]
    raise ValueError(f"Invalid input type: {config['input_type']}")


def get_mapping_inform(read_key: str, mode: str) -> Path:
    """
    Returns the mapping informs.
    :param read_key: Read key
    :param mode: Reference mode ('assembly' or 'ref')
    :return: Path to mapping informs
    """
    if read_key == 'fastq_pe':
        if mode == 'assembly':
            return Path(str(OUTPUT_ASSEMBLY_MAPPING_INFORMS).format(mapper='bowtie2'))
        elif mode == 'ref':
            return Path(str(OUTPUT_VARIANT_CALLING_MAPPING_INFORMS))
    elif read_key == 'fastq_se':
        if mode == 'ref':
            logger.warning(f'Reference mapping for single-end data is not implemented yet')
        return Path(str(OUTPUT_ASSEMBLY_MAPPING_INFORMS).format(mapper='minimap2'))
    else:
        raise ValueError(f'Invalid read key: {read_key}')


def get_depth_inform(read_key: str, mode: str) -> Path:
    """
    Returns the depth informs.
    :param read_key: Read key
    :param mode: Reference mode ('assembly' or 'ref')
    :return: Path to depth informs
    """
    if read_key == 'fastq_pe':
        if mode == 'assembly':
            return Path(str(OUTPUT_ASSEMBLY_DEPTH_INFORMS).format(mapper='bowtie2'))
        else:
            return Path(str(OUTPUT_VARIANT_CALLING_DEPTH_INFORMS))
    elif read_key == 'fastq_se':
        if mode == 'ref':
            logger.warning(f'Reference mapping for single-end data is not implemented yet')
        return Path(str(OUTPUT_ASSEMBLY_DEPTH_INFORMS).format(mapper='minimap2'))
    else:
        raise ValueError(f'Invalid read key: {read_key}')


def get_mapping_rate_inform(read_key: str) -> Path:
    """
    Returns the mapping rate informs.
    :param read_key: Read key
    :return: Path to depth informs
    """
    if read_key == 'fastq_pe':
        return Path(str(OUTPUT_ASSEMBLY_MAPPING_RATE_INFORMS).format(mapper='bowtie2'))
    elif read_key == 'fastq_se':
        return Path(str(OUTPUT_ASSEMBLY_MAPPING_RATE_INFORMS).format(mapper='minimap2'))
    else:
        raise ValueError(f'Invalid read key: {read_key}')


def get_qc_informs(config: Dict[str, Any], input_type: str, mode: str = 'assembly') -> List[Path]:
    """
    Returns the QC informs based on the input type.
    :param config: Snakemake configuration
    :param input_type: Input type
    :param mode: Reference mode (assembly or ref)
    :return: List of paths to QC informs
    """
    informs = []
    if input_type in ('hybrid', 'illumina'):
        informs.append(get_mapping_inform('fastq_pe', mode))
        informs.append(get_depth_inform('fastq_pe', mode))
    if input_type in ('hybrid', 'ont'):
        informs.append(get_mapping_inform('fastq_se', mode))
        informs.append(get_depth_inform('fastq_se', mode))
    return [Path(config['working_dir'], p) for p in informs]
