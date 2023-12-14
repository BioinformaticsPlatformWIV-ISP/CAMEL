import unittest
from pathlib import Path
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.tools.pipelines.salmonella.spifinder import SPIFinder
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliodirectory import ToolIODirectory


class TestSPIFinderFasta(CamelTestSuite):
    """
    Tests the SPIFinder tool for fasta input
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('salmonella')
    input_fasta_file = test_file_dir / 'assembly_filtered.fasta'

    def test_spifinder_fasta(self) -> None:
        """
        Tests basic spifinder run on fasta input data.
        :return: None
        """
        spifinder = SPIFinder(self.camel)
        spifinder.add_input_files({
            'FASTA': [ToolIOFile(Path(TestSPIFinderFasta.input_fasta_file))],
            'DIR': [ToolIODirectory(Path('/db/pipelines/salmonella/spifinder/genomicepidemiology-'
                                         'spifinder_db-db102668b704'))]
        })
        spifinder.run(self.running_dir)
        self.verify_output_files(spifinder, 'JSON')


if __name__ == '__main__':
    unittest.main()
