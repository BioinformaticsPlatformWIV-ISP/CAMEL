import copy
import shutil
import tempfile
from pathlib import Path

from Bio.Phylo import NewickIO
from Bio.Phylo.Newick import Tree

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.figtree.figtree import FigTree
from camel.app.tools.treevector.treevector import TreeVector


class NewickUtils(object):
    """
    This class contains utility functions to work with phylogenetic trees in Newick format.
    """

    @staticmethod
    @PendingDeprecationWarning
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
        logger.info(f"Rendering tree: {tree_path}")
        tree_vector = TreeVector(camel)
        tree_vector.add_input_files({'NWK': [ToolIOFile(tree_path)]})
        tree_vector.update_parameters(output_format=output_format, output_filename='tree.png', type=plot_type)
        tree_vector.run(Path(tempfile.mkdtemp(prefix='tree_vector', dir=camel.config['temp_dir'])))
        logger.debug(f"TreeVector - stdout: {tree_vector.stdout}")
        logger.debug(f"TreeVector - stderr: {tree_vector.stderr}")
        shutil.copyfile(tree_vector.tool_outputs['PNG'][0].path, output_path)

    @staticmethod
    def calculate_tree_image_height(min_height: int, nb_leaves: int) -> int:
        """
        Calculates the required height for the given tree.
        :param min_height: Minimum height
        :param nb_leaves: Nb. of leaves in the tree
        :return: Image height
        """
        return min_height + 20 * nb_leaves

    @staticmethod
    def create_image_figtree(path_newick: Path, path_config: Path, path_out: Path, width: int, height: int) -> FigTree:
        """
        Creates an image for the given Newick tree using FigTree.
        :param path_newick: Path to Newick file
        :param path_config: Path to FigTree config file
        :param path_out: Output path
        :param width: Image width
        :param height: Image height
        :return: FigTree instance
        """
        logger.info(f"Creating image for tree: {path_newick}")
        with tempfile.TemporaryDirectory(prefix='figtree_', dir=Camel.get_instance().config['temp_dir']) as dir_temp:
            figtree = FigTree(Camel.get_instance())
            output_path = Path(dir_temp, 'tree.png')
            figtree.update_parameters(output_path=str(path_out), width=width, height=height)
            path_template = path_config
            figtree.add_input_files({'NWK': [ToolIOFile(path_newick)], 'TXT': [ToolIOFile(path_template)]})
            figtree.run(Path(dir_temp))
        logger.info(f"PNG visualization created: {output_path}")
        return figtree

    @staticmethod
    def count_leaves(newick_in: Path) -> int:
        """
        Counts the number of leaves in a phylogenetic tree.
        :param newick_in: Input Newick file
        :return: Number of leaves
        """
        with open(newick_in) as handle:
            tree = next(NewickIO.parse(handle))
        return tree.count_terminals()

    @staticmethod
    def remove_inner_node_names(tree: Tree) -> Tree:
        """
        Removes inner node names from a tree.
        :param tree: Tree object
        :return: Updated tree
        """
        logger.info("Removing inner node names")
        new_tree = copy.copy(tree)
        for clade in new_tree.find_clades():
            if 'Inner' in clade.name:
                clade.name = ''
        return new_tree

    @staticmethod
    def export_newick_tree(input_tree: Tree, path_out: Path) -> None:
        """
        Saves a tree to a file in Newick format.
        :param input_tree: Input tree
        :param path_out: Output path
        :return: None
        """
        logger.info(f"Saving tree to file: {path_out.name}")
        with path_out.open('w') as handle:
            NewickIO.write([input_tree], handle)
