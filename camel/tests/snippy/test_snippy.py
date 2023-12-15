import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.snippy.snippy import Snippy


class TestSnippy(CamelTestSuite):
    """
    Tests the snippy tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('snippy')
    input_fasta_ref = test_file_dir / 'H37Rv.fasta'
    input_fasta_ctgs = test_file_dir / 'S32BD03301_contigs_unfilt.fasta'

    def test_snippy(self) -> None:
        """
        Tests the snippy tool
        :return: None
        """
        snippy = Snippy(self.camel)
        snippy.add_input_files({'FASTA_CTGS': [ToolIOFile(Path(TestSnippy.input_fasta_ctgs))],
                                'FASTA_REF': [ToolIOFile(Path(TestSnippy.input_fasta_ref))]
                                })
        snippy.run(self.running_dir)
        self.verify_output_files(snippy, 'BAM', 1)


if __name__ == '__main__':
    unittest.main()