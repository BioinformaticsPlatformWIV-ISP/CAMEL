from pathlib import Path
from typing import Dict, Any

from camel.resources.snakefile import assembly_spades, assembly_flye


def get_fasta(config: Dict[str, Any]) -> Path:
    """
    Returns the assembly FASTA output IO object path.
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
