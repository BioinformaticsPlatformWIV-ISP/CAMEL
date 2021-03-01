#!/usr/bin/env python
import argparse
import logging
from typing import Dict, List, Tuple, Sequence, Optional

from Bio.Phylo.Newick import Tree
from Bio.Phylo.TreeConstruction import DistanceMatrix

from camel.app.camel import Camel
from camel.app.components.phylogeny.mlstphylotutils import MlstPyhloUtils
from camel.app.components.phylogeny.mlstreportparser import MlstReportParser
from camel.app.components.phylogeny.mlsttabularparser import MlstTabularParser
from camel.app.components.phylogeny.newickutils import NewickUtils


class MainMlstTree(object):
    """
    The main class for the MLST tree tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main scripts.
        :param args: Command line arguments
        """
        self._args = MainMlstTree._parse_arguments(args)

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: (optional) arguments
        :return: Arguments
        """
        ap = argparse.ArgumentParser()
        ap.add_argument('--input-html', nargs=2, action='append')
        ap.add_argument('--input-tab', nargs=2, action='append')
        ap.add_argument('--clustering-method', choices=['nj', 'upgma'], default='nj')
        ap.add_argument('--output', type=str)
        ap.add_argument('--output-image', type=str)
        ap.add_argument('--output-tabular', type=str)
        ap.add_argument('--output-dist-matrix', type=str)
        ap.add_argument('--plot-type', default='clad', choices=['clad', 'phylo'])
        ap.add_argument('--include-imperfect-hits', action='store_true')
        ap.add_argument('--no-tree', action='store_true',
                        help="If set, no output tree is generated. Useful when the input data is really big.")
        return ap.parse_args(args)

    def run(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        allele_ids_by_sample = self.__parse_input_files()

        # Create the tabular output file
        if self._args.output_tabular:
            self.__create_tabular_output(allele_ids_by_sample, self._args.output_tabular)

        # Stop when only tabular output is specified
        if self._args.output is None and self._args.output_dist_matrix is None:
            logging.info("Not creating tree (output not set)")
            return

        # Construct distance matrix
        matrix = MlstPyhloUtils.calculate_distance_matrix(allele_ids_by_sample)
        if self._args.output_dist_matrix is not None:
            self.__export_distance_matrix(matrix)

        # Create tree
        if self._args.no_tree is not True:
            tree = MlstPyhloUtils.construct_tree(matrix, self._args.clustering_method)
            self.__export_tree(tree)

    def __create_tabular_output(self, allele_ids_by_sample: Dict[str, List[Tuple[str, str]]], output_path: str) -> None:
        """
        Creates a tabular output file with the detected alleles for all of the samples.
        :param allele_ids_by_sample: Allele ids by sample
        :param output_path: Output path
        :return: None
        """
        locus_names = self.__get_locus_names(allele_ids_by_sample)

        # Create tabular output
        sample_names = sorted(list(allele_ids_by_sample.keys()))
        header = ['Locus'] + sample_names
        table_data = [header]
        for i in range(0, len(locus_names)):
            row = [locus_names[i]]
            for sample_name in sample_names:
                row.append(allele_ids_by_sample[sample_name][i][1])
            table_data.append(row)

        # Save to file
        logging.info(f"Creating tabular output file: {output_path}")
        with open(output_path, 'w') as handle:
            for row in table_data:
                handle.write('\t'.join(row))
                handle.write('\n')

    def __get_locus_names(self, allele_ids_by_sample: Dict[str, List[Tuple[str, str]]]) -> List[str]:
        """
        Returns the locus names.
        A loop is used to return only the first element as the locus names are identical for all entries in the
        dictionary. An error is raised when the allele id dictionary is empty.
        :param allele_ids_by_sample: Allele ids by sample
        :return: List of locus names
        """
        for _, allele_ids_by_sample in allele_ids_by_sample.items():
            return [locus for locus, _ in allele_ids_by_sample]
        raise ValueError("Cannot determine locus names")

    def __parse_input_files(self) -> Dict[str, List[Tuple[str, str]]]:
        """
        Parses the input files for the tree construction.
        :return: Detected alleles by sample name
        """
        allele_ids_by_sample = {}
        if self._args.input_html:
            allele_ids_by_sample = MlstReportParser.parse_html_all(
                self._args.input_html, self._args.include_imperfect_hits)
        elif self._args.input_tab:
            allele_ids_by_sample = MlstTabularParser.parse_tabular_all(
                self._args.input_tab, self._args.include_imperfect_hits)
        if len(allele_ids_by_sample) < 3:
            raise ValueError("At least 3 samples are required")
        return allele_ids_by_sample

    def __export_distance_matrix(self, matrix: DistanceMatrix) -> None:
        """
        Exports the distance matrix to a text file.
        :param matrix: Distance matrix
        :return: None
        """
        with open(self._args.output_dist_matrix, 'w') as h_out:
            h_out.write(str(matrix))
        logging.info(f"Distance matrix exported to: {self._args.output_dist_matrix}")

    def __export_tree(self, tree: Tree) -> None:
        """
        Exports the generated tree.
        :return: None
        """
        camel = Camel()
        tree = NewickUtils.remove_inner_node_names(tree)
        NewickUtils.export_newick_tree(tree, self._args.output)
        if self._args.output_image:
            NewickUtils.render(camel, self._args.output, self._args.output_image, self._args.plot_type)


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    main = MainMlstTree()
    main.run()
