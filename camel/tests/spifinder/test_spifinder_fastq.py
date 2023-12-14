import unittest
from pathlib import Path
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.tools.pipelines.salmonella.spifinder import SPIFinder
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliodirectory import ToolIODirectory


class TestSPIFinderFastq(CamelTestSuite):
    """
    Tests the SPIFinder tool for fastq input
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('salmonella')
    input_pe_reads = [test_file_dir / "SRR493330_1.fastq.gz", test_file_dir / "SRR493330_2.fastq.gz"]

    def test_spifinder_fastq(self) -> None:
        """
        Tests basic spifinder run on fastq input data.
        :return: None
        """
        spifinder = SPIFinder(self.camel)
        spifinder.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in self.input_pe_reads],
            'DIR': [ToolIODirectory(Path('/db/pipelines/salmonella/spifinder/genomicepidemiology-'
                                         'spifinder_db-db102668b704'))]
        })
        spifinder.run(self.running_dir)
        self.verify_output_files(spifinder, 'JSON')


if __name__ == '__main__':
    unittest.main()
