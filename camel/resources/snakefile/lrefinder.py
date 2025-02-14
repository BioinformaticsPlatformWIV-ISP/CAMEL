from pathlib import Path
from typing import Any

from camel.resources.snakefile import read_simulation

SNAKEFILE_LREFINDER = f'{Path(__file__).parent / Path(__file__).stem}.smk'
OUTPUT_LREFINDER_REPORT = Path('lrefinder') / 'html.io'
OUTPUT_LREFINDER_REPORT_EMPTY = Path('lrefinder') / 'html_empty.io'
OUTPUT_LREFINDER_INFORMS = Path('lrefinder') / 'informs.io'
OUTPUT_LREFINDER_SUMMARY = Path('lrefinder') / 'summary.tsv'


def get_input(config: dict[str, Any]) -> Path:
    """
    Returns the input for the LREFinder tool.
    :param config: Snakemake configuration
    :return: Path to the input file
    """
    if config['input_type'] in ('illumina', 'ont', 'hybrid'):
        return Path(config['working_dir']) / 'fq_dict.io'
    if config['input_type'] == 'fasta':
        return Path(config['working_dir']) / read_simulation.OUTPUT_SIMULATION_FASTQ
    raise ValueError(f"Input type '{config['input_type']}' is not supported by LRE-finder")