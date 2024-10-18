import json
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

SNAKEFILE_SPIFINDER = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_spifinder = Path('spifinder')
OUTPUT_SPIFINDER_FASTQ_JSON = _dir_spifinder / 'spifinder_fastq' / '{input_format}' / 'spifinder_output.io'
OUTPUT_SPIFINDER_FASTA_JSON = _dir_spifinder / 'spifinder_fasta' / 'spifinder_output.io'
OUTPUT_SPIFINDER_REPORT = _dir_spifinder / '{input_format}' / 'html.io'
OUTPUT_SPIFINDER_FASTQ_INFORMS = _dir_spifinder / 'spifinder_fastq' / '{input_format}' / 'informs.io'
OUTPUT_SPIFINDER_FASTA_INFORMS = _dir_spifinder / 'spifinder_fasta' / 'informs.io'
OUTPUT_SPIFINDER_REPORT_EMPTY = _dir_spifinder / '{input_format}' / 'html-empty.io'
OUTPUT_SPIFINDER_SUMMARY = _dir_spifinder / '{input_format}' / 'summary_out.tsv'
OUTPUT_SPIFINDER_DOC = _dir_spifinder / '{input_format}' / 'spifinder_function_category.tsv'
OUTPUT_SPIFINDER_SUMMARY_IO = _dir_spifinder / '{input_format}' / 'summary_out.io'


def get_reports(config: Dict[str, Any]) -> List[Path]:
    """
    Returns the paths to the spifinder reports.
    :param config: Snakemake configuration
    :return: Report path(s)
    """
    input_type = config['input_type']
    paths = []

    if input_type == 'fasta':
        if 'spifinder' in config['analyses']:
            paths.append(str(OUTPUT_SPIFINDER_REPORT).format(input_format='fasta'))
        else:
            paths.append(str(OUTPUT_SPIFINDER_REPORT_EMPTY).format(input_format='fasta'))

    # PE reads
    if input_type in ('illumina', 'hybrid'):
        if 'spifinder' in config['analyses']:
            paths.append(str(OUTPUT_SPIFINDER_REPORT).format(input_format='fastq_pe'))
        else:
            paths.append(str(OUTPUT_SPIFINDER_REPORT_EMPTY).format(input_format='fastq_pe'))

    # SE reads
    if input_type in ('ont', 'hybrid'):
        if 'spifinder' in config['analyses']:
            paths.append(str(OUTPUT_SPIFINDER_REPORT).format(input_format='fastq_se'))
        else:
            paths.append(str(OUTPUT_SPIFINDER_REPORT_EMPTY).format(input_format='fastq_se'))

    return [Path(config['working_dir']) / p for p in paths]


def get_command_informs(config: Dict[str, Any]) -> List[Path]:
    """
    Returns a list of informs IO files for spifinder.
    :param config: Snakemake configuration
    :return: List of informs IO files
    """
    input_type = config['input_type']
    paths = []

    if 'spifinder' not in config['analyses']:
        return []

    # FASTA input
    if input_type == 'fasta':
        paths.append(OUTPUT_SPIFINDER_FASTA_INFORMS)

    # PE reads
    if input_type in ('illumina', 'hybrid'):
        paths.append(OUTPUT_SPIFINDER_FASTA_INFORMS)
        paths.append(str(OUTPUT_SPIFINDER_FASTQ_INFORMS).format(input_format='fastq_pe'))

    # SE reads
    if input_type in ('ont', 'hybrid'):
        paths.append(OUTPUT_SPIFINDER_FASTA_INFORMS)
        paths.append(str(OUTPUT_SPIFINDER_FASTQ_INFORMS).format(input_format='fastq_se'))
    return [Path(config['working_dir']) / p for p in paths]


def get_summaries(config: Dict[str, Any]) -> List[Path]:
    """
    Returns the paths to the spifinder summary file(s).
    :param config: Snakemake configuration
    :return: Summary file path(s)
    """
    input_type = config['input_type']
    paths = []

    if 'spifinder' not in config['analyses']:
        return []

    # FASTA input
    if input_type == 'fasta':
        paths.append(str(OUTPUT_SPIFINDER_SUMMARY).format(input_format='fasta'))

    # PE reads
    if input_type in ('illumina', 'hybrid'):
        paths.append(str(OUTPUT_SPIFINDER_SUMMARY).format(input_format='fastq_pe'))

    # SE reads
    if input_type in ('ont', 'hybrid'):
        paths.append(str(OUTPUT_SPIFINDER_SUMMARY).format(input_format='fastq_se'))

    return [Path(config['working_dir']) / p for p in paths]


def spifinder_json_parser(json_file_path: Path, mode: str) -> List[Union[str, int, float]]:
    """
    This function is able to parse the output json files of the spifinder tool and returns more favorable outputs
    for the TSV.
    :param json_file_path: Path of the json file to be parsed
    :param mode: fasta or fastq
    :return: a list of hits to be added in the output tsv, a dictionary of hits and metadata to be added to the
    output json
    """
    with json_file_path.open('r') as file_handle:
        json_file = json.load(file_handle)

    # define header for TSV:
    header_part1 = ['SPI', 'identity', 'coverage']
    if mode == 'fasta':
        header_part2 = ['contig_name', 'positions_in_contig', 'accession', 'insertion_site',
                        'category_function']
    elif mode == 'fastq':
        header_part2 = ['accession', 'insertion_site', 'category_function']
    else:
        raise ValueError(f"This function's parameter 'mode' must be either fastq or fasta, current value is {mode}")

    results_tsv = []
    spi = json_file['spifinder']["results"]['Salmonella Pathogenicity Islands']['SPI']
    if not spi == "No hit found":
        for hit in spi.keys():
            # Add coverage to spi hit dictionary because it doesn't exist
            spi[hit]['coverage'] = f"{spi[hit]['HSP_length']}/{spi[hit]['template_length']}"

            for hit_property in header_part1 + header_part2:
                # add to tsv
                results_tsv.append(spi[hit][hit_property])

    return results_tsv
