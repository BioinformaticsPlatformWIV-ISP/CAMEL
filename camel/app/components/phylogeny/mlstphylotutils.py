import logging
from typing import List, Dict, Tuple

from Bio.Phylo.BaseTree import Tree
# noinspection PyProtectedMember
from Bio.Phylo.TreeConstruction import _DistanceMatrix, DistanceTreeConstructor


class MlstPyhloUtils(object):
    """
    This class contains utility functions to generate MLST based phylogenetic trees.
    """

    @staticmethod
    def calculate_distance(alleles_a: List[Tuple[str, str]], alleles_b: List[Tuple[str, str]]) -> float:
        """
        Calculates the distance between two sets of detected alleles.
        :param alleles_a: Gene names + allele ids of first sample
        :param alleles_b: Gene names + allele ids of second sample
        :return: Calculated distance between samples
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

            if allele_id_a not in ('-', '?') and allele_id_b not in ('-', '?'):
                genes_common.append(gene_name_a)

                if not allele_id_a == allele_id_b:
                    genes_different.append(gene_name_a)
        return float(len(genes_different)) / len(genes_common) if len(genes_common) > 0 else 1

    @staticmethod
    def calculate_distance_matrix(allele_ids: Dict[str, List[Tuple[str, str]]]) -> _DistanceMatrix:
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
                row.append(MlstPyhloUtils.calculate_distance(allele_ids[sample_names[i]], allele_ids[sample_names[j]]))
            row.append(0)
            matrix.append(row)
        return _DistanceMatrix(sample_names, matrix)

    @staticmethod
    def construct_tree(matrix: _DistanceMatrix, clustering_method: str) -> Tree:
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
