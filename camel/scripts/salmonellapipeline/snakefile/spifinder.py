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
    results = []
    spi = json_file['spifinder']["results"]['Salmonella Pathogenicity Islands']['SPI']
    if spi == "No hit found":
        inter_json_dict = {f"spifinder_{mode}": {'results': results}}
    else:
        hit_dictionary_list = []
        for hits in spi.keys():
            json_dict = {}
            header_part1 = ['SPI', 'identity']
            if mode == 'fasta':
                header_part2 = ['contig_name', 'positions_in_contig', 'accession', 'insertion_site',
                                'category_function']
            else:  # mode == 'fastq':
                header_part2 = ['accession', 'insertion_site', 'category_function']
            results.append([spi[hits][hit_property] for hit_property in header_part1] +
                           [f"{spi[hits]['HSP_length']}/{spi[hits]['template_length']}"] +
                           [spi[hits][hit_property] for hit_property in header_part2])
            for hit_property in header_part1 + header_part2:
                json_dict[hit_property] = spi[hits][hit_property]
            json_dict['coverage'] = f"{spi[hits]['HSP_length']}/{spi[hits]['template_length']}"
            hit_dictionary_list.append(json_dict)
        inter_json_dict = {f"spifinder_{mode}": {'results': hit_dictionary_list}}

    inter_json_dict[f"spifinder_{mode}"]['informs_tools'] = {
        tool_informs.get('_tool', tool_informs['_name']): {'_name': tool_informs['_name'],
                                                           '_version': tool_informs['_version'],
                                                           '_command': tool_informs['_command'],
                                                           '_tag': tool_informs['_tag']}}
    inter_json_dict[f"spifinder_{mode}"]['informs_dbs'] = {'last_updated': tool_informs['last_update_date'],
                                                           'name': tool_informs['key'], 'title': tool_informs['key']}
    return results, inter_json_dict
