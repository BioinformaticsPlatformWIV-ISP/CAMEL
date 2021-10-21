import copy
import logging
import os
import shutil
import tempfile
from pathlib import Path

from Bio.Phylo import NewickIO
from Bio.Phylo.Newick import Tree

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.treevector.treevector import TreeVector


class NewickUtils(object):
    """
    This class contains utility functions to work with phylogenetic trees in Newick format.
    """

    @staticmethod
    def render(camel: Camel, tree_path: Path, output_path: Path, plot_type: str, output_format: str = 'png') -> None:
        """
        Renders the tree to image.
        :param camel: Camel instance
        :param tree_path: Tree file in Newick format
        :param output_path: Image output
        :param plot_type: Type of plot
        :param output_format: Image output format
        :return: None
        """
        logging.info(f"Rendering tree: {tree_path}")
        tree_vector = TreeVector(camel)
        tree_vector.add_input_files({'NWK': [ToolIOFile(tree_path)]})
        tree_vector.update_parameters(output_format=output_format, output_filename='tree.png', type=plot_type)
        tree_vector.run(Path(tempfile.mkdtemp(prefix='tree_vector', dir=camel.config['temp_dir'])))
        logging.debug(f"TreeVector - stdout: {tree_vector.stdout}")
        logging.debug(f"TreeVector - stderr: {tree_vector.stderr}")
        shutil.copyfile(tree_vector.tool_outputs['PNG'][0].path, output_path)

    @staticmethod
    def remove_inner_node_names(tree: Tree) -> Tree:
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

    @staticmethod
    def export_newick_tree(input_tree: Tree, filename: str) -> None:
        """
        Saves a tree to a file in Newick format.
        :param input_tree: Input tree
        :param filename: Filename
        :return: None
        """
        logging.info(f"Saving tree to file: {os.path.basename(filename)}")
        with open(filename, 'w') as handle:
            NewickIO.write([input_tree], handle)
