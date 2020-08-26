from camel.scripts.influenzapipeline.snakefile import subtyping_blastn
from typing import Dict, Any, Union, List
from pathlib import Path


def get_subtyping_report(config: Dict[str, Any]) -> Union[Path, List[Any]]:
    """
    Returns the path to the subtyping report.
    :param config: Snakemake configuration
    :return: Path to subtyping report
    """
    if config['subtyping_method'] == 'blastn':
        return Path(config['working_dir']) / subtyping_blastn.OUTPUT_SUBTYPING_REPORT
    else:
        return []


def get_subtyping_summary(config: Dict[str, Any]) -> Union[Path, List[Any]]:
    """
    Returns the path to the subtyping summary.
    :param config: Snakemake configuration
    :return: Path to subtyping summary
    """
    if config['subtyping_method'] == 'blastn':
        return Path(config['working_dir']) / subtyping_blastn.OUTPUT_SUBTYPING_SUMMARY
    else:
        return []
