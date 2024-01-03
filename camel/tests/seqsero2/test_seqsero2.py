import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.pipelines.salmonella.seqsero2 import SeqSero2


class TestSeqsero2(CamelTestSuite):
    """
    Tests the Seqsero2 tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('salmonella')
    input_pe_reads = [test_file_dir / "SRR493330_1.fastq.gz", test_file_dir / "SRR493330_2.fastq.gz"]
    input_fasta_file = test_file_dir / 'assembly_filtered.fasta'

    def test_seqsero2_kmer(self) -> None:
        """
        Tests basic seqsero2 run in Kmer mode.
        :return: None
        """
        seqserotool = SeqSero2(self.camel)
        seqserotool.add_input_files({
            'FASTA': [ToolIOFile(Path(self.input_fasta_file))],
            'DIR': [ToolIODirectory(Path('/db/pipelines/salmonella/seqsero2/1.2.1/seqsero2_db'))],
            'MODE': [ToolIOValue('Kmer')]
        })
        seqserotool.run(self.running_dir)
        self.verify_output_files(seqserotool, 'TXT')

    def test_seqsero2_allele(self) -> None:
        """
        Tests basic seqsero2 run in Allele mode.
        :return: None
        """
        seqserotool = SeqSero2(self.camel)
        seqserotool.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in self.input_pe_reads],
            'DIR': [ToolIODirectory(Path('/db/pipelines/salmonella/seqsero2/1.2.1/seqsero2_db'))],
            'MODE': [ToolIOValue('Allele')]
        })
        seqserotool.run(self.running_dir)
        self.verify_output_files(seqserotool, 'TXT')

    def test_seqsero2_kmerread(self) -> None:
        """
        Tests basic seqsero2 run in Kmerread mode.
        :return: None
        """
        seqserotool = SeqSero2(self.camel)
        seqserotool.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in self.input_pe_reads],
            'DIR': [ToolIODirectory(Path('/db/pipelines/salmon  ella/seqsero2/1.2.1/seqsero2_db'))],
            'MODE': [ToolIOValue('Kmerread')]
        })
        seqserotool.run(self.running_dir)
        self.verify_output_files(seqserotool, 'TXT')


if __name__ == '__main__':
    unittest.main()
