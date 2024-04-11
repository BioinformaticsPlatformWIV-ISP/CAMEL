import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.art.art import ART


class TestART(CamelTestSuite):
    """
    Tests the ART tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('snippy')
    input_fasta = test_file_dir / 'S32BD03301_contigs_unfilt.fasta'

    def test_art(self) -> None:
        """
        Tests the ART tool.
        :return: None
        """
        art = ART(self.camel)
        art.add_input_files({
            'FASTA': [ToolIOFile(Path(TestART.input_fasta))]
        })
        art.run(self.running_dir)
        self.verify_output_files(art, 'FASTQ_PE', 2)


if __name__ == '__main__':
    unittest.main()
