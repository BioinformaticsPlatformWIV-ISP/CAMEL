import argparse
import logging
from typing import Dict, List, Tuple

from camel.app.camel import Camel
from camel.scripts.mlsttree import htmlinputparser, treeconstruction, tabularinputparser


def _parse_arguments():
    """
    Parses the command line arguments.
    :return: Arguments
    """
    ap = argparse.ArgumentParser()
    ap.add_argument('--input-html', nargs=2, action='append')
    ap.add_argument('--input-tab', nargs=2, action='append')
    ap.add_argument('--clustering-method', choices=['nj', 'upgma'])
    ap.add_argument('--output')
    ap.add_argument('--output-image')
    ap.add_argument('--output-tabular')
    return ap.parse_args()


def __create_tabular_output(a_ids: Dict[str, List[Tuple[str, str]]], output_path: str) -> None:
    """
    Creates a tabular output file with the detected alleles for all of the samples.
    :param a_ids: Allele ids by sample
    :param output_path: Output path
    :return: None
    """
    locus_names = __get_locus_names(a_ids)

    # Create tabular output
    sample_names = sorted(list(a_ids.keys()))
    header = ['Locus'] + sample_names
    table_data = [header]
    for i in range(0, len(locus_names)):
        row = [locus_names[i]]
        for sample_name in sample_names:
            row.append(a_ids[sample_name][i][1])
        table_data.append(row)

    # Save to file
    logging.info(f"Creating tabular output file: {output_path}")
    with open(args.output_path, 'w') as handle:
        for row in table_data:
            handle.write('\t'.join(row))
            handle.write('\n')


def __get_locus_names(a_ids: Dict[str, List[Tuple[str, str]]]) -> List[str]:
    """
    Returns the locus names.
    :param a_ids: Allele ids by sample
    :return: List of locus names
    """
    for _, a_ids in a_ids.items():
        return [locus for locus, _ in a_ids]
    raise ValueError("Cannot determine locus names")


if __name__ == '__main__':
    camel = Camel()
    args = _parse_arguments()
    allele_ids_by_sample = {}
    if args.input_html:
        allele_ids_by_sample = htmlinputparser.parse_all(args.input_html)
    elif args.input_tab:
        allele_ids_by_sample = tabularinputparser.parse_all(args.input_tab)
    if len(allele_ids_by_sample) < 3:
        raise ValueError("At least 3 samples are required")

    if args.output_tabular:
        __create_tabular_output(allele_ids_by_sample, args.output_tabular)

    matrix = treeconstruction.calculate_distance_matrix(allele_ids_by_sample)
    tree = treeconstruction.generate(matrix, args.clustering_method)
    tree = treeconstruction.remove_inner_node_names(tree)
    treeconstruction.save_tree(tree, args.output)
    if args.output_image:
        treeconstruction.render(camel, args.output, args.output_image)
