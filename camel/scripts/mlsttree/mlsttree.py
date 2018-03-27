import argparse

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
    return ap.parse_args()


if __name__ == '__main__':
    camel = Camel()
    args = _parse_arguments()
    allele_ids = {}
    if args.input_html:
        allele_ids = htmlinputparser.parse_all(args.input_html)
    elif args.input_tab:
        allele_ids = tabularinputparser.parse_all(args.input_tab)
    if len(allele_ids) < 3:
        raise ValueError("At least 3 samples are required")

    matrix = treeconstruction.calculate_distance_matrix(allele_ids)
    tree = treeconstruction.generate(matrix, args.clustering_method)
    tree = treeconstruction.remove_inner_node_names(tree)
    treeconstruction.save_tree(tree, args.output)
    if args.output_image:
        treeconstruction.render(camel, args.output, args.output_image)
