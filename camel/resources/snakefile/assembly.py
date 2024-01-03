from pathlib import Path
from typing import Dict, Any

from camel.resources.snakefile import assembly_spades, assembly_flye

SNAKEFILE_ASSEMBLY = f'{Path(__file__).parent / Path(__file__).stem}.smk'
OUTPUT_ASSEMBLY_FASTA = Path('assembly', 'filtering', 'fasta.io')
OUTPUT_ASSEMBLY_FILTERING_INFORMS = Path('assembly', 'filtering', 'informs.io')
OUTPUT_ASSEMBLY_REPORT = Path('assembly', 'report', 'html.io')
OUTPUT_ASSEMBLY_SUMMARY = Path('assembly', 'summary', 'summary_out.tsv')
OUTPUT_ASSEMBLY_MAPPING_INFORMS = Path('assembly', '{mapper}', 'informs.io')
OUTPUT_ASSEMBLY_DEPTH_INFORMS = Path('assembly', 'samtools_depth', '{mapper}', 'informs.io')


def get_fasta_raw(config: Dict[str, Any]) -> Path:
    """
    Returns the assembly FASTA output IO object path (before filtering).
    """
    if config['input_type'] == 'illumina':
        return assembly_spades.OUTPUT_ASSEMBLY_FASTA
    if config['input_type'] in ('ont', 'hybrid'):
        return assembly_flye.OUTPUT_ASSEMBLY_FASTA
    raise ValueError(f"Invalid input type: {config['input_type']}")


def get_command_informs(config: Dict[str, Any]) -> Path:
    """
    Returns the assembly informs output IO object path.
    :return: Assembly informs path
    """
    if config['input_type'] == 'illumina':
        return Path(config['working_dir'], assembly_spades.OUTPUT_ASSEMBLY_INFORMS)
    if config['input_type'] in ('ont', 'hybrid'):
        return Path(config['working_dir'], assembly_flye.OUTPUT_ASSEMBLY_INFORMS)
    raise ValueError(f"Invalid input type: {config['input_type']}")


def get_mapping_inform(read_key: str) -> Path:
    """
    Returns the mapping informs.
    :param read_key: Read key
    :return: Path to mapping informs
    """
    if read_key == 'fastq_pe':
        return Path(str(OUTPUT_ASSEMBLY_MAPPING_INFORMS).format(mapper='bowtie2'))
    elif read_key == 'fastq_se':
        return Path(str(OUTPUT_ASSEMBLY_MAPPING_INFORMS).format(mapper='minimap2'))
    else:
        raise ValueError(f'Invalid read key: {read_key}')


def get_depth_inform(read_key: str) -> Path:
    """
    Returns the depth informs.
    :param read_key: Read key
    :return: Path to depth informs
    """
    if read_key == 'fastq_pe':
        return Path(str(OUTPUT_ASSEMBLY_DEPTH_INFORMS).format(mapper='bowtie2'))
    elif read_key == 'fastq_se':
        return Path(str(OUTPUT_ASSEMBLY_DEPTH_INFORMS).format(mapper='minimap2'))
    else:
        raise ValueError(f'Invalid read key: {read_key}')
