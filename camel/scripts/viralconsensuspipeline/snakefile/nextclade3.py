from pathlib import Path
from typing import Any

from camel.app.snakemake.snakemakeutils import SnakemakeUtils

SNAKEFILE_NEXTCLADE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_nextclade = Path('nextclade')
OUTPUT_NEXTCLADE_SUBTYPE_REPORT = _dir_nextclade / 'subtype_determination' / 'report' / 'html.io'
OUTPUT_NEXTCLADE_SUBTYPE_REPORT_EMPTY = _dir_nextclade / 'subtype_determination' / 'report' / 'html-empty.io'
OUTPUT_NEXTCLADE_REPORT = _dir_nextclade / 'html.io'
OUTPUT_NEXTCLADE_REPORT_EMPTY = _dir_nextclade / 'html-empty.io'
OUTPUT_NEXTCLADE_SUMMARY = _dir_nextclade / 'summary_nextclade.tsv'
OUTPUT_NEXTCLADE_INFORMS = _dir_nextclade / 'informs.io'


def get_nextclade_db(wildcards, cps, segment: str, config: dict) -> str:
    """
    Returns the path to the Nextclade database.
    :param wildcards: Snakemake wildcards
    :param cps: Snakemake checkpoints
    :param segment: Segment
    :param config: Snakemake config
    :return: Path to the database
    """
    if config['nextclade'].get('dbs') is not None:
        return config['nextclade']['dbs'][segment]

    # Get data from subtype determination
    # noinspection PyUnresolvedReferences
    path_informs = Path(cps.nextclade3_detect_subtype_report.get().output['INFORMS'])
    informs_subtype = SnakemakeUtils.load_object(path_informs)
    return informs_subtype['nextclade_dbs'][segment]


def get_informs_subtype(wildcards: Any, cps: Any) -> Path:
    """
    Returns the subtype determination informs.
    :param wildcards: Snakemake wildcards
    :param cps: Snakemake checkpoints
    :return: Path to subtype informs
    """
    path_informs = Path(cps.nextclade3_detect_subtype_report.get().output['INFORMS'])
    return path_informs


def get_nextclade_output(wildcards, cps: Any, key: str, config: dict = None) -> list[Path]:
    """
    Aggregates the Nextclade output based on the database information.
    :param wildcards: Snakemake wildcards
    :param cps: Snakemake checkpoints
    :param config: Configuration
    :param key: Output key (TSV / INFORMS)
    :return: List of Nextclade outputs
    """
    # Extract segments
    if config['nextclade'].get('dbs') is not None:
        segments = list(config['nextclade']['dbs'].keys())
    else:
        # noinspection PyUnresolvedReferences
        path_informs = Path(cps.nextclade3_detect_subtype_report.get().output['INFORMS'])
        informs_subtype = SnakemakeUtils.load_object(path_informs)
        segments = list(informs_subtype['nextclade_dbs'].keys())

    # Determine output
    if key == 'TSV':
        base_output = Path('nextclade', '{segment}', 'tsv.io')
    elif key == 'INFORMS':
        base_output = Path('nextclade', '{segment}', 'informs.io')
    else:
        raise ValueError(f'Invalid key: {key}')

    # Return list of outputs
    return [Path(config['working_dir'], str(base_output).format(segment=segment)) for segment in segments]
