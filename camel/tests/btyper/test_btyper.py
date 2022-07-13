import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.btyper.btyper import BTyper

class TestBtyper(CamelTestSuite):


    """
    Initializes this testing tool
    """
    test_file_dir = Path('/testdata/camel/btyper/')
    FILE_FASTA = ToolIOFile(test_file_dir / 'bacillus_contigs.fasta')

    def test_btyper(self) -> None:
        """
        actually testing BTyper
        """
        btyper = BTyper(self.camel)
        btyper.add_input_files({'FASTA':[TestBtyper.FILE_FASTA]})
        btyper.update_parameters(output_dir = self.running_dir)
        btyper.run(self.running_dir)
        self.verify_output_files(btyper, 'TSV')


if __name__ == '__main__':
    unittest.main()