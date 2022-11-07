import os
from pathlib import Path
from typing import Dict, List

from snakemake.checkpoints import Checkpoints
from snakemake.io import expand

SNAKEFILE_FLEXNERITYPE = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_flexneri = Path('flexneri_type')
OUTPUT_FLEXNERI_REPORT = _dir_flexneri / 'report' / 'val-html.io'
OUTPUT_FLEXNERI_REPORT_EMPTY = _dir_flexneri / 'report' / 'val-html-empty.io'
OUTPUT_FLEXNERI_SUMMARY = _dir_flexneri / 'summary_out.tsv'


def aggregate_input(wildcards: Dict, checkpoints: Checkpoints, config: Dict) -> List[str]:
    """
    Aggregates the input from the read mapping based on the detected loci.
    :param wildcards: Wildcards
    :param checkpoints: Checkpoints
    :param config: Snakemake configuration
    :return: List of input files
    """
    dir_fasta = checkpoints.flexneri_type_prepare_reference_files.get(**wildcards).output.DIR_FASTA
    loci = [locus for locus in os.listdir(dir_fasta) if not locus.startswith('.')]
    return expand(str(Path(config['working_dir']) / 'flexneri_type' / 'loci' / '{flexneri_locus}' / 'val-mut.io'),
                  flexneri_locus=loci)


def aggregate_input_vcf(wildcards: Dict, checkpoints: Checkpoints, config: Dict) -> List[Path]:
    """
    Aggregates the input from the read mapping based on the detected loci.
    :param wildcards: Wildcards
    :param checkpoints Checkpoints
    :param config: Snakemake configuration
    :return: List of input files
    """
    dir_fasta = checkpoints.flexneri_type_prepare_reference_files.get(**wildcards).output.DIR_FASTA
    loci = [locus for locus in os.listdir(dir_fasta) if not locus.startswith('.')]
    return expand(str(Path(config['working_dir']) / 'flexneri_type' / 'loci' / '{flexneri_locus}' / 'vcf-csq.io'),
                  flexneri_locus=loci)
