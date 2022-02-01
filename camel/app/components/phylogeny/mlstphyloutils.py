import json
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional

import pandas as pd


def is_perfect(record: pd.Series, detection_method: str) -> bool:
    """
    Determines if the given hit is perfect.
    :param record: Input record
    :param detection_method: Allele detection method.
    :return: True if perfect, False otherwise
    """
    if detection_method == 'blast':
        if str(record['Allele']) in ('-', '?'):
            return False
        if float(record['% Identity']) != 100.0:
            return False
        len_hsp, len_locus = record['HSP/Locus length'].split('/')
        if not len_hsp == len_locus:
            return False
        return True
    elif detection_method == 'srst2':
        if str(record['Allele']) == '-':
            return False
        if str(record['Mismatches']) != '0':
            return False
        if str(record['Uncertainty']) != '-':
            return False
        return True
    elif detection_method == 'kma':
        if str(record['Allele']) == '-':
            return False
        if float(record['% Identity']) != 100.0:
            return False
        if float(record['% Coverage']) != 100.0:
            return False
        return True
    else:
        raise ValueError(f"Invalid detection method: {detection_method}")


def parse_tsv_typing(tsv_path: Path, detection_method: str) -> Dict[str, str]:
    """
    Parses a tabular output file for the sequence typing assay.
    :param tsv_path: Typing output file
    :param detection_method: Allele detection method
    :return: Parsed alleles (key: locus, value: allele as a string)
    """
    allele_data = pd.read_table(tsv_path)
    allele_data['is_perfect_hit'] = allele_data.apply(lambda x: is_perfect(x, detection_method), axis=1)
    return {r['Locus']: r['Allele'] if r['is_perfect_hit'] else '-' for _, r in allele_data.iterrows()}


def parse_tsv_typing_list(tsv_in: List[Tuple[Path, str]], detection_method: Optional[str] = 'blast') -> pd.DataFrame:
    """
    Parses a list of tabular input files.
    :param tsv_in: List of input files + file name
    :param detection_method: Detection method for sequence typing
    :return: Dictionary of detected alleles by sample name
    """
    sample_names = []
    allele_data = []
    for tabular_file, file_name in tsv_in:
        logging.debug(f'Parsing file: {tabular_file}')
        try:
            isolate_name = Path(file_name).stem
            logging.debug(f'Sample name: {isolate_name}')
        except IndexError:
            raise ValueError(f'Cannot determine sample name from: {file_name}')
        allele_data.append(parse_tsv_typing(tabular_file, detection_method))
        sample_names.append(isolate_name)
    return pd.DataFrame(allele_data, index=sample_names, dtype=str)


def parse_html_typing_list(dirs_in: List[Path], detection_method: Optional[str] = 'blast') -> pd.DataFrame:
    """
    Parses a list of HTML output directories.
    :param dirs_in: List of input directories
    :param detection_method: Detection method for sequence typing
    :return: Dictionary of detected alleles by sample name
    """
    sample_names = []
    allele_data = []
    for dir_typing in dirs_in:
        # Retrieve metadata from JSON file
        json_typing_metadata = list((dir_typing / 'sequence_typing').glob('**/sequence_typing.json'))
        if len(json_typing_metadata) == 0:
            raise FileNotFoundError("Cannot find typing metadata")
        elif len(json_typing_metadata) != 1:
            raise FileNotFoundError("Multiple typing outputs found")
        with json_typing_metadata[0].open() as handle:
            typing_metadata = json.load(handle)

        # Add sample data to allele data
        sample_names.append(typing_metadata['sample'])
        tsv_typing = next(json_typing_metadata[0].parent.glob('typing-*.tsv'))
        allele_data.append(parse_tsv_typing(tsv_typing, detection_method))
    return pd.DataFrame(allele_data, index=sample_names, dtype=str)
