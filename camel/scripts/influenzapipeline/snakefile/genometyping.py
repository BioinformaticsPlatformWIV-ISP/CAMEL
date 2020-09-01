from pathlib import Path
from typing import Dict, Any, Union, List

from camel.scripts.influenzapipeline.snakefile import genometyping_blastn


def get_segmenttyping_report(config: Dict[str, Any]) -> Union[Path, List[Any]]:
    """
    Returns the path to the genometyping report.
    :param config: Snakemake configuration
    :return: Path to genometyping report
    """
    if config['genometyping_method'] == 'blastn':
        return Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_REPORT
    else:
        return []


def get_genometyping_summary(config: Dict[str, Any]) -> Union[Path, List[Any]]:
    """
    Returns the path to the genometyping summary.
    :param config: Snakemake configuration
    :return: Path to genometyping summary
    """
    if config['genometyping_method'] == 'blastn':
        return Path(config['working_dir']) / genometyping_blastn.OUTPUT_GENOMETYPING_SUMMARY
    else:
        return []
