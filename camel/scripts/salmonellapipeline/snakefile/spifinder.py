import json
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

SNAKEFILE_SPIFINDER = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_spifinder = Path('spifinder')
OUTPUT_JSON_SPIFINDER_FASTQ = _dir_spifinder / 'spifinder_fastq' / 'spifinder_output.io'
OUTPUT_JSON_SPIFINDER_FASTA = _dir_spifinder / 'spifinder_fasta' / 'spifinder_output.io'
OUTPUT_SPIFINDER_REPORT = _dir_spifinder / 'html.io'
OUTPUT_SPIFINDER_FASTQ_INFORMS = _dir_spifinder / 'spifinder_fastq' / 'informs.io'
OUTPUT_SPIFINDER_FASTA_INFORMS = _dir_spifinder / 'spifinder_fasta' / 'informs.io'
OUTPUT_SPIFINDER_REPORT_EMPTY = _dir_spifinder / 'html-empty.io'
OUTPUT_SPIFINDER_SUMMARY = _dir_spifinder / 'summary_out.tsv'
OUTPUT_SPIFINDER_SUMMARY_JSON = _dir_spifinder / 'summary_out.json'
OUTPUT_SPIFINDER_DOC = _dir_spifinder / 'spifinder_function_category.tsv'
OUTPUT_SPIFINDER_SUMMARY_IO = _dir_spifinder / 'summary_out.io'


def spifinder_json_parser(json_file_path: Path, tool_informs: Dict[str, Any], mode: str) -> \
        Tuple[List[List[Union[str, int, float], ...]], Dict[Any]]:
    """
    This function is able to parse the output json files of the spifinder tool and returns more favorable outputs
    for the Tsv and the camel Json for Hera.
    :param json_file_path: Path of the json file to be parsed
    :param tool_informs: tool informs corresponding to the run of which json_file_path was the output
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
    if spi == "No hit found":
        inter_json_dict = {f"spifinder_{mode}": {'results': results_tsv}}
    else:
        hit_dictionary_list = []
        for hit in spi.keys():
            json_dict = {}

            # Add coverage to spi hit dictionary because it doesn't exist
            spi[hit]['coverage'] = f"{spi[hit]['HSP_length']}/{spi[hit]['template_length']}"

            for hit_property in header_part1 + header_part2:
                # add to tsv
                results_tsv.append(spi[hit][hit_property])
                # add to json
                json_dict[hit_property] = spi[hit][hit_property]

            hit_dictionary_list.append(json_dict)
        inter_json_dict = {f"spifinder_{mode}": {'results': hit_dictionary_list}}

    # Add tool and db metadata to json dictionary
    inter_json_dict[f"spifinder_{mode}"]['informs_tools'] = {
        tool_informs.get('_tool', tool_informs['_name']): {'_name': tool_informs['_name'],
                                                           '_version': tool_informs['_version'],
                                                           '_command': tool_informs['_command'],
                                                           '_tag': tool_informs['_tag']}}
    inter_json_dict[f"spifinder_{mode}"]['informs_dbs'] = {'last_updated': tool_informs['last_update_date'],
                                                           'name': tool_informs['key'], 'title': tool_informs['key']}
    return results_tsv, inter_json_dict
