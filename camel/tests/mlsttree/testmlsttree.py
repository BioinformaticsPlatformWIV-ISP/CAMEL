import itertools
import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.scripts.mlsttree.mainmlsttree import MainMlstTree


class TestMlstTree(CamelTestSuite):
    """
    Tests the MLST tree tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('mlst_tree')
    input_tabular_files = [
        test_file_dir / 'typing-cgmlst-S1.tsv',
        test_file_dir / 'typing-cgmlst-S2.tsv',
        test_file_dir / 'typing-cgmlst-S3.tsv',
        test_file_dir / 'typing-cgmlst-S4.tsv',
        test_file_dir / 'typing-cgmlst-S5.tsv'
    ]

    def test_tree_construction(self) -> None:
        """
        Tests the tree construction.
        :return: None
        """
        output_file_newick = self.running_dir / 'my_tree.nwk'
        output_file_tabular = self.running_dir / 'my_tree.tsv'
        output_file_dist_matrix = self.running_dir / 'dist_matrix.txt'
        output_file_image = self.running_dir / 'my_tree.png'
        args = [
            '--output-image', str(output_file_image),
            '--output-tabular', str(output_file_tabular),
            '--output-dist-matrix', str(output_file_dist_matrix),
            '--output', str(output_file_newick)
        ] + list(
            itertools.chain.from_iterable([['--input-tab', str(p), p.name] for p in TestMlstTree.input_tabular_files]))
        mlst_tree = MainMlstTree(args)
        mlst_tree.run()
        self.assertGreater(output_file_newick.stat().st_size, 0)
        self.assertGreater(output_file_tabular.stat().st_size, 0)
        self.assertGreater(output_file_dist_matrix.stat().st_size, 0)
        self.assertGreater(output_file_image.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
