import logging

import os


def __parse_tabular_output(output_file, file_name):
    """
    Parses a tabular output file.
    :param output_file: Output file
    :param file_name: Output file name
    :return: Allele ids
    """
    with open(output_file) as handle:
        header = handle.readline()
        if not header.split('\t')[0] == 'Locus' and header.split('\t')[1] == 'Allele':
            raise ValueError("Invalid tabular file: {}".format(file_name))
        allele_ids = [(line.split('\t')[0], line.split('\t')[1]) for line in handle.readlines()]
        logging.info('Nb. of alleles: {}'.format(len(allele_ids)))
        return allele_ids


def parse(tabular_input):
    """
    Parses the tabular input.
    :param tabular_input: Tabular input
    :return: Allele ids
    """
    allele_ids = {}
    for tabular_file, file_name in tabular_input:
        logging.info('Parsing file: {}'.format(file_name))
        try:
            sample_name = os.path.splitext(file_name)[0]
            logging.info('Sample: {}'.format(sample_name))
        except IndexError:
            raise ValueError("Cannot determine sample name from: {}".format(file_name))
        allele_ids[sample_name] = __parse_tabular_output(tabular_file, sample_name)
    return allele_ids
