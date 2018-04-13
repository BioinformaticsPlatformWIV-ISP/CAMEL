import copy
import logging

import shutil
from Bio.Phylo import NewickIO

# noinspection PyProtectedMember
from Bio.Phylo.TreeConstruction import _DistanceMatrix, DistanceTreeConstructor

from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.treevector.treevector import TreeVector


def __calculate_distance(alleles_a, alleles_b):
    """
    Calculates the distance between two samples.
    :param alleles_a: First sample
    :param alleles_b: Second sample
    :return: Distance
    """
    if not len(alleles_a) == len(alleles_b):
        raise ValueError('Length of input data does not match.')
    genes_common = []
    genes_different = []

    for i in range(0, len(alleles_a)):
        gene_name_a, allele_id_a = alleles_a[i]
        gene_name_b, allele_id_b = alleles_b[i]

        if not gene_name_a == gene_name_b:
            raise ValueError(f'Gene names do not match ({gene_name_a}, {gene_name_b})')

        if not allele_id_a == '-' and not allele_id_b == '-':
            genes_common.append(gene_name_a)

            if not allele_id_a == allele_id_b:
                genes_different.append(gene_name_a)
    return float(len(genes_different)) / len(genes_common) if len(genes_common) > 0 else 1


def calculate_distance_matrix(allele_ids):
    """
    Calculates a distance matrix.
    :param allele_ids: The allele ids in the format {sample_name: [[gene_name, allele], ...]}.
    :return: distance matrix, sample names
    """
    logging.info("Calculating distance matrix")
    sample_names = sorted(allele_ids.keys())
    matrix = [[0]]
    for j in range(1, len(sample_names)):
        row = []
        for i in range(0, j):
            row.append(__calculate_distance(allele_ids[sample_names[i]], allele_ids[sample_names[j]]))
        row.append(0)
        matrix.append(row)
    return _DistanceMatrix(sample_names, matrix)


def generate(matrix, clustering_method):
    """
    Generates a tree for the given distance matrix using the given clustering method.
    :param matrix: Distance matrix
    :param clustering_method: Clustering method
    :return: Tree
    """
    logging.info("Generating tree")
    constructor = DistanceTreeConstructor()
    if clustering_method == 'upgma':
        return constructor.upgma(matrix)
    elif clustering_method == 'nj':
        return constructor.nj(matrix)
    else:
        raise ValueError('Unknown clustering method "{}"'.format(clustering_method))


def remove_inner_node_names(tree):
    """
    Removes inner node names from a tree
    :param tree: Tree
    :return:
    """
    logging.info("Removing inner node names")
    new_tree = copy.copy(tree)
    for clade in new_tree.find_clades():
        if 'Inner' in clade.name:
            clade.name = ''
    return new_tree


def save_tree(input_tree, filename):
    """
    Saves a tree to a file.
    :param input_tree: Input tree
    :param filename: Filename
    :return: None
    """
    logging.info("Saving tree to file")
    with open(filename, 'w') as handle:
        NewickIO.write([input_tree], handle)


def render(camel, tree_file, output):
    """
    Renders the tree to image.
    :param camel: Camel instance
    :param tree_file: Tree file in Newick format
    :param output: Image output
    :return: None
    """
    logging.info("Rendering tree")
    tree_vector = TreeVector(camel)
    tree_vector.add_input_files({'NWK': [ToolIOFile(tree_file)]})
    tree_vector.update_parameters(output_format='png', output_filename='tree.png', type='clad')
    tree_vector.run()
    shutil.copyfile(tree_vector.tool_outputs['PNG'][0].path, output)
