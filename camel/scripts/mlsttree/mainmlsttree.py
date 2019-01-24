#!/usr/bin/env python
import argparse
import logging
from typing import Dict, List, Tuple

from Bio.Phylo.Newick import Tree

from camel.app.camel import Camel
from camel.app.components.phylogeny.mlsphylotutils import MlstPyhloUtils
from camel.app.components.phylogeny.mlstreportparser import MlstReportParser
from camel.app.components.phylogeny.mlsttabularparser import MlstTabularParser
from camel.app.components.phylogeny.newickutils import NewickUtils


class MainMlstTree(object):
    """
    The main class for the MLST tree tool.
    """

    def __init__(self, args: argparse.Namespace=None):
        """
        Initializes the main scripts.
        :param args: Command line arguments
        """
        self._args = args if args is not None else MainMlstTree._parse_arguments()

    @staticmethod
    def _parse_arguments() -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Arguments
        """
        ap = argparse.ArgumentParser()
        ap.add_argument('--input-html', nargs=2, action='append')
        ap.add_argument('--input-tab', nargs=2, action='append')
        ap.add_argument('--clustering-method', choices=['nj', 'upgma'])
        ap.add_argument('--output', type=str)
        ap.add_argument('--output-image', type=str)
        ap.add_argument('--output-tabular', type=str)
        ap.add_argument('--plot-type', default='clad', choices=['clad', 'phylo'])
        return ap.parse_args()

    def run(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        allele_ids_by_sample = self.__parse_input_files()
        if self._args.output_tabular:
            self.__create_tabular_output(allele_ids_by_sample, self._args.output_tabular)
        matrix = MlstPyhloUtils.calculate_distance_matrix(allele_ids_by_sample)
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
            allele_ids_by_sample = MlstReportParser.parse_html_all(self._args.input_html)
        elif self._args.input_tab:
            allele_ids_by_sample = MlstTabularParser.parse_tabular_all(self._args.input_tab)
        if len(allele_ids_by_sample) < 3:
            raise ValueError("At least 3 samples are required")
        return allele_ids_by_sample

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
