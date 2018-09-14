import argparse
import unittest

import os
import tempfile

from camel.app.camel import Camel
from camel.scripts.mlsttree.mainmlsttree import MainMlstTree


class TestMlstTree(unittest.TestCase):
    """
    Tests the MLST tree tool.
    """
    camel = Camel()
    running_dir = None

    # Input files
    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'])
    input_tabular_files = [
        os.path.join(test_file_dir, 'mlst_tree', 'typing-cgmlst-S1.tsv'),
        os.path.join(test_file_dir, 'mlst_tree', 'typing-cgmlst-S2.tsv'),
        os.path.join(test_file_dir, 'mlst_tree', 'typing-cgmlst-S3.tsv'),
        os.path.join(test_file_dir, 'mlst_tree', 'typing-cgmlst-S4.tsv'),
        os.path.join(test_file_dir, 'mlst_tree', 'typing-cgmlst-S5.tsv')
    ]

    def setUp(self) -> None:
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(prefix='camel_', dir=TestMlstTree.camel.config['temp_dir'])

    def test_tree_construction(self) -> None:
        """
        Tests the tree construction.
        :return: None
        """
        output_file_newick = os.path.join(self.running_dir, 'my_tree.nwk')
        output_file_tabular = os.path.join(self.running_dir, 'my_tree.tsv')
        output_file_image = os.path.join(self.running_dir, 'my_tree.png')
        args = argparse.Namespace(
            input_tab=[(p, os.path.basename(p)) for p in TestMlstTree.input_tabular_files],
            input_html=None,
            output_image=output_file_image,
            output_tabular=output_file_tabular,
            output=output_file_newick,
            clustering_method='upgma',
            plot_type='clad',
            working_dir=self.running_dir
        )
        mlst_tree = MainMlstTree(args)
        mlst_tree.run()
        self.assertGreater(os.path.getsize(output_file_newick), 0)
        self.assertGreater(os.path.getsize(output_file_tabular), 0)
        self.assertGreater(os.path.getsize(output_file_image), 0)
