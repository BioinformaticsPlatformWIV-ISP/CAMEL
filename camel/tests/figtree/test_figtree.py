import unittest
from pathlib import Path

import pkg_resources

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.figtree.figtree import FigTree


class TestFigTree(CamelTestSuite):
    """
    Tests the FigTree tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('newick')
    input_nwk = test_file_dir / 'tree_mst.nwk'

    def test_figtree(self) -> None:
        """
        Tests the FigTree tool.
        :return: None
        """
        path_png_out = self.running_dir / 'tree.png'
        figtree = FigTree(self.camel)
        figtree.add_input_files({'NWK': [ToolIOFile(TestFigTree.input_nwk)]})
        figtree.update_parameters(output_path=str(path_png_out))
        figtree.run()
        self.assertIn('PNG', figtree.tool_outputs)
        self.assertGreater(figtree.tool_outputs['PNG'][0].size, 0)

    def test_figtree_with_template(self) -> None:
        """
        Tests the FigTree tool with a template.
        :return: None
        """
        path_png_out = self.running_dir / 'tree.png'
        path_template = Path(pkg_resources.resource_filename('camel', 'resources/figtree/template_cgmlst_tree.txt'))
        figtree = FigTree(self.camel)
        figtree.add_input_files({'NWK': [ToolIOFile(TestFigTree.input_nwk)], 'TXT': [ToolIOFile(path_template)]})
        figtree.update_parameters(output_path=str(path_png_out))
        figtree.run()
        self.assertIn('PNG', figtree.tool_outputs)
        self.assertGreater(figtree.tool_outputs['PNG'][0].size, 0)


if __name__ == '__main__':
    unittest.main()
