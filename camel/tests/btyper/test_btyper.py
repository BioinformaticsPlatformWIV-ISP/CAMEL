import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.btyper.btyper import Btyper

class TestBtyper(CamelTestSuite):
    """
    initializes this testing tool
    """

    test_file_dir = Path('/home/magodfroid/PycharmProjects/CamelProjects/tmp')
    FILE_FASTA = ToolIOFile(test_file_dir / 'test_contigs.fasta')

    def test_btyper(self) -> None:
        """
        actually testing btyper
        """

        btyper = Btyper(self.camel)
        btyper.add_input_files({'FASTA':[TestBtyper.FILE_FASTA]})
        btyper.update_parameters(output_dir='tmp2')
        btyper.run(self.running_dir)

if __name__ == '__main__':
    unittest.main()