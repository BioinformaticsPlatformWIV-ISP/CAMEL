import logging

import json
import os


def __get_scheme_dir(html_file, html_name):
    """
    Returns the directory with the output of the MLST scheme.
    :param html_file: HTML file
    :param html_name: HTML Galaxy name
    :return: Scheme directory.
    """
    files_dir = f"{html_file.split('.')[0]}_files"
    st_dir = os.path.join(files_dir, 'sequence_typing')
    if not os.path.isdir(st_dir):
        raise ValueError("Input file '{}' is not a valid sequence typing output".format(html_name))
    if len(os.listdir(st_dir)) != 1:
        raise ValueError("Multiple sequence typing outputs found")
    return os.path.join(st_dir, os.listdir(st_dir)[0])


def __get_analysis_info(scheme_dir):
    """
    Returns the analysis metadata from the given scheme directory.
    :return: Analysis info dictionary
    """
    info_path = os.path.join(scheme_dir, 'analysis_info.txt')
    if not os.path.isfile(info_path):
        raise ValueError("No analysis info file found (you might have to rerun the sequence typing tool)")
    with open(info_path) as handle:
        metadata = json.load(handle)
        logging.info('Scheme: {}'.format(metadata['scheme']))
        logging.info('Sample: {}'.format(metadata['sample']))
    return metadata


def __parse_st_output(scheme_dir):
    """
    Parses the sequence typing output file.
    :return: Allele ids
    """
    try:
        tabular_file = [os.path.join(scheme_dir, x) for x in os.listdir(scheme_dir) if x.endswith('.tsv')][0]
    except IndexError:
        raise FileNotFoundError("No tabular output file found")
    with open(tabular_file) as handle:
        header = handle.readline()
        if header.split('\t')[0] != 'Locus' or header.split('\t')[1] != 'Allele':
            raise ValueError("Invalid format in tabular file: {}".format(tabular_file))
        allele_ids = [(line.split('\t')[0], line.split('\t')[1]) for line in handle.readlines()]
        logging.info('Nb. of alleles: {}'.format(len(allele_ids)))
        return allele_ids


def parse_all(html_input):
    """
    Parses all HTML input files that were provided trough the command line arguments.
    :return: Allele ids
    """
    allele_ids = {}
    for html_file, html_name in html_input:
        logging.info('Parsing file: {}'.format(html_name))
        scheme_dir = __get_scheme_dir(html_file, html_name)
        metadata = __get_analysis_info(scheme_dir)
        if metadata['sample'] in allele_ids:
            logging.warning("Duplicate sample! {}".format(metadata['sample']))
            continue
        allele_ids[metadata['sample']] = __parse_st_output(scheme_dir)
        logging.debug('-'*10)
    return allele_ids
