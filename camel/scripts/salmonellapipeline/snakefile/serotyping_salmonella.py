import re
from pathlib import Path
from typing import Any, Dict, List, Tuple

SNAKEFILE_SEROTYPE = f'{Path(__file__).parent / Path(__file__).stem}.smk'
_dir_serotype = Path('serotyping')
OUTPUT_SISTR_SEROTYPE = _dir_serotype / 'serotyping_sistr' / 'sistr_output.io'
OUTPUT_SEROTYPE_SISTR_INFORMS = _dir_serotype / 'serotyping_sistr' / 'informs.io'
OUTPUT_SEQSERO2_SEROTYPE_KMER_INFORMS = _dir_serotype / 'serotyping_seqsero2_Kmer' / 'informs.io'
OUTPUT_SEQSERO2_SEROTYPE_ALLELE_INFORMS = _dir_serotype / 'serotyping_seqsero2_Allele' / 'informs.io'
OUTPUT_SEQSERO2_SEROTYPE_KMERREAD_INFORMS = _dir_serotype / 'serotyping_seqsero2_Kmerread' / 'informs.io'
OUTPUT_SEROTYPE_REPORT = _dir_serotype / 'html.io'
OUTPUT_SEROTYPE_REPORT_EMPTY = _dir_serotype / 'html-empty.io'
OUTPUT_SEROTYPE_SUMMARY = _dir_serotype / 'summary_out.tsv'
OUTPUT_SEROTYPE_SUMMARY_JSON = _dir_serotype / 'summary_out.json'


def sistr_output_parser(prediction: Dict[str, Any], locus: str, antigen: str, hits_dict_tsv: Dict[str, str], 
                        hits_dict_json: Dict[str, Any], header_locus: List[str]) -> None:
    """
    Parses the Sistr output for a specific locus (the o antigen has 2 loci, and the h antigens each have one locus).
    Updates the hits dictionaries in place without returning any output.
    :param prediction: the dictionary of the results of the specific locus
    :param locus: locus name, either fliC, fljB, wzx, or wzy
    :param antigen: antigen name, either h1, h2, or o
    :param hits_dict_tsv: the dictionary of the results for the tsv file
    :param hits_dict_json: the nested dictionary of the results for the json file
    :param header_locus: the header for both the tsv and json dictionary values
    :return: None
    """
    is_missing = str(prediction['is_missing'])
    hit_properties = [
        locus,
        prediction[antigen if antigen != 'o' else 'serogroup'].replace(',', ';'),
        format(prediction['top_result']['pident'], '.2f'),
        '/'.join([str(prediction['top_result']['length']),
                  str(prediction['top_result']['qlen'])]),
        prediction['top_result']['stitle'],
        '...'.join([str(prediction['top_result']['sstart']),
                    str(prediction['top_result']['send'])])
    ]

    if is_missing:
        json_dict = {item: "-" for item in header_locus}
    else:
        json_dict = {header_locus[i]: hit_properties[i] for i in range(len(header_locus))}

    hits_dict_tsv[f'hits_serotype_{antigen}_{locus}'] = ','.join(hit_properties)
    hits_dict_json[f'hits_serotype_{antigen}_{locus}'] = json_dict


def seqsero2_output_parser(seqsero2_file: Path, seqsero2_mode: str, informs_dict: Dict[str, str]) -> \
        Tuple[Dict[str, Any], List[str]]:
    """
    Parses the output file of a SeqSero2 run, uniform over all three modes.
    :param seqsero2_file: path of the output file
    :param seqsero2_mode: mode of the SeqSero2 run, either seqsero2_kmer, seqsero2_allele, or seqsero2_kmerread
    :param informs_dict: corresponding informs dictionary of the given seqsero2_file and mode
    :return: tuple of 1. intermediate dictionary to be combined and then written to json file and
    2. List of result string to be written to tsv file
    """
    json_dict = {}
    with seqsero2_file.open('r') as handle:
        tsv_results = handle.readlines()[2:8]
    tsv_results = [re.sub(r'([^ ]) ([^ ])', r'\1_\2', res).strip("\n") for res in tsv_results]
    tsv_results = [res.replace(':\t', '\t') for res in tsv_results]
    tsv_results = [seqsero2_mode + "_" + x for x in tsv_results]
    for res in tsv_results:
        json_dict[res.split('\t')[0]] = res.split('\t')[1]
    inter_json_dict = {
        seqsero2_mode: {
            **json_dict,
            'informs_tools': {informs_dict.get('_tool', informs_dict['_name']): {
                '_name': informs_dict['_name'], '_version': informs_dict['_version'],
                '_command': informs_dict['_command'], '_tag': informs_dict['_tag']}},
            'informs_dbs': {
                'last_updated': informs_dict['last_update_date'], 'name': informs_dict['key'],
                'title': informs_dict['key']}}}
    return inter_json_dict, tsv_results
