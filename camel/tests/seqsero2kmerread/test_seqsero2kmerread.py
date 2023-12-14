import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.pipelines.salmonella.seqsero2kmerread import SeqSero2KmerRead


class TestSeqsero2kmerread(CamelTestSuite):
    """
    Tests the Seqsero2 tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('salmonella')
    input_pe_reads = [test_file_dir / "SRR493330_1.fastq.gz", test_file_dir / "SRR493330_2.fastq.gz"]

    def test_seqsero2kmerread(self) -> None:
        """
        Tests basic seqsero2 run.
        :return: None
        """

        seqserotool = SeqSero2KmerRead(self.camel)
        seqserotool.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in self.input_pe_reads],
            'DIR': [ToolIODirectory(Path('/db/pipelines/salmonella/seqsero2/1.2.1/seqsero2_db'))]
        })
        seqserotool.run(self.running_dir)
        self.verify_output_files(seqserotool, 'TXT')


if __name__ == '__main__':
    unittest.main()
