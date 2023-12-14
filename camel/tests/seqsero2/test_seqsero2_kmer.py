import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.pipelines.salmonella.seqsero2 import SeqSero2


class TestSeqsero2Kmer(CamelTestSuite):
    """
    Tests the Seqsero2 tool in
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('salmonella')
    input_fasta_file = test_file_dir / 'assembly_filtered.fasta'

    def test_seqsero2_kmer(self) -> None:
        """
        Tests basic seqsero2 run.
        :return: None
        """
        seqserotool = SeqSero2(self.camel)
        seqserotool.add_input_files({
            'FASTA': [ToolIOFile(Path(TestSeqsero2Kmer.input_fasta_file))],
            'DIR': [ToolIODirectory(Path('/db/pipelines/salmonella/seqsero2/1.2.1/seqsero2_db'))]
        })
        seqserotool.run(self.running_dir)
        self.verify_output_files(seqserotool, 'TXT')


if __name__ == '__main__':
    unittest.main()
