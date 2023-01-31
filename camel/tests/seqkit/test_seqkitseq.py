import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.seqkit.seqkitseq import SeqkitSeq


class TestSeqkitSeq(CamelTestSuite):
    """
    Initializes the SeqKit seq testing tool
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('seqkit')
    FILE_FASTA = ToolIOFile(test_file_dir / 'assembly-VAR305.fasta')
    FILE_FASTQ = ToolIOFile(test_file_dir / 'reads_illumina_1.fastq')

    def test_seqkitseq_fasta(self) -> None:
        """
        Testing SeqKit seq with contigs file.
        """
        seqkitseq = SeqkitSeq(self.camel)
        seqkitseq.add_input_files({'FASTA': [TestSeqkitSeq.FILE_FASTA]})
        seqkitseq.update_parameters(output_filename='sequences.fasta', max_length=100000)
        seqkitseq.run(self.running_dir)
        self.verify_output_files(seqkitseq, 'FASTA')

    def test_seqkitseq_fastq(self) -> None:
        """
        Testing SeqKit seq with paired-end fastq reads.
        """
        seqkitseq = SeqkitSeq(self.camel)
        seqkitseq.add_input_files({'FASTQ': [TestSeqkitSeq.FILE_FASTQ]})
        seqkitseq.update_parameters(max_length=150)
        seqkitseq.run(self.running_dir)
        self.verify_output_files(seqkitseq, 'FASTQ')


if __name__ == '__main__':
    unittest.main()
