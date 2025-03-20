from pathlib import Path
from typing import Any, Dict, List

SNAKEFILE_STRAINGST = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_straingst = Path('straingst', '{read_type}')
OUTPUT_HDF5_STRAINGST = _dir_straingst / 'hdf5-straingst.io'
OUTPUT_INFORMS_STRAINGST = _dir_straingst / 'informs.io'
OUTPUT_STRAINGST_REPORT = _dir_straingst / 'html.io'
OUTPUT_STRAINGST_REPORT_EMPTY = _dir_straingst / 'html-empty.io'
OUTPUT_STRAINGST_SUMMARY = _dir_straingst / 'summary_out.tsv'
OUTPUT_STRAINGST_STATS = _dir_straingst / 'straingst-stats.io'
OUTPUT_STRAINGST_STRAINS = _dir_straingst / 'straings-strains.io'


def get_command_informs(config: Dict[str, Any]) -> List[Path]:
    """
    Returns a list of informs IO files for the straingst steps.
    :param config: Snakemake configuration
    :return: List of informs IO files
    """
    paths = []

    # Straingst is disabled -> return empty list
    if not any([an for an in config['analyses'] if an in ['straingst', 'gmo']]):
        return []

    input_type = config['input_type']
    if input_type in ('illumina', 'hybrid'):
        paths.append(Path(str(OUTPUT_INFORMS_STRAINGST).format(read_type='illumina')))
    if input_type == 'ont':
        paths.append(Path(str(OUTPUT_INFORMS_STRAINGST).format(read_type='ont')))
    return [Path(config['working_dir']) / p for p in paths]


def get_reports(config: Dict[str, Any]) -> List[Path]:
    """
    Returns the path to the strainGST report.
    :param config: Snakemake configuration
    :return: List of paths to the straingst reports
    """
    paths = []
    input_type = config['input_type']

    # FASTA input
    if config['input_type'] == 'fasta':
        return []

    # FASTQ input
    if input_type in ('illumina', 'hybrid'):
        if any([an for an in config['analyses'] if an in ['straingst', 'gmo']]):
            paths.append(Path(str(OUTPUT_STRAINGST_REPORT).format(read_type='illumina')))
        else:
            paths.append(Path(str(OUTPUT_STRAINGST_REPORT_EMPTY).format(read_type='illumina')))
    if input_type == 'ont':
        if any([an for an in config['analyses'] if an in ['straingst', 'gmo']]):
            paths.append(Path(str(OUTPUT_STRAINGST_REPORT).format(read_type='ont')))
        else:
            paths.append(Path(str(OUTPUT_STRAINGST_REPORT_EMPTY).format(read_type='ont')))

    # Check if at least one path was added
    if len(paths) == 0:
        raise ValueError(f'Invalid input type: {input_type}')

    return [Path(config['working_dir']) / p for p in paths]


def get_summaries(config: Dict[str, Any]) -> List[Path]:
    """
    Returns the paths to the straingst summary file(s).
    :param config: Snakemake configuration
    :return: List of paths to straingst summary files
    """
    input_type = config['input_type']

    # FASTA input
    if config['input_type'] == 'fasta':
        return []

    # FASTQ input
    paths = []
    # StrainGST is for now ran on Illumina reads due to parallelization issues
    # Also, can be part of the GMO assay or can also be ran independently, hence the following conditions cover that
    # scenario.
    if (input_type in ('illumina', 'hybrid')) and any([an for an in config['analyses'] if an in ['straingst', 'gmo']]):
        paths.append(Path(str(OUTPUT_STRAINGST_SUMMARY).format(read_type='illumina')))
    if (input_type == 'ont') and any([an for an in config['analyses'] if an in ['straingst', 'gmo']]):
        paths.append(Path(str(OUTPUT_STRAINGST_SUMMARY).format(read_type='ont')))
    return [Path(config['working_dir']) / p for p in paths]
