import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.scripts.sequencetyping.mainsequencetyping import MainSequenceTyping


class TestSequenceTyping(CamelTestSuite):
    """
    Tests the sequence typing tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('typing')
    input_db_nucl = test_file_dir / 'scheme_mlst_neisseria'
    input_db_protein = test_file_dir / 'scheme_pora_neisseria'
    input_db_mixed = test_file_dir / 'scheme_fhbp_neisseria'
    input_fasta = test_file_dir / 'neisseria_mc58.fasta'
    input_typing_reads = [
        test_file_dir / 'S15BD05018_S58_L001_1.fastq',
        test_file_dir / 'S15BD05018_S58_L001_2.fastq'
    ]

    def test_typing_blast_nucl(self) -> None:
        """
        Tests sequence typing using BLAST with a nucleotide scheme (including ST definitions).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta),
            '--scheme-dir', str(self.input_db_nucl),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_blast_pept(self) -> None:
        """
        Tests sequence typing using BLAST with a peptide scheme.
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta),
            '--scheme-dir', str(self.input_db_protein),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8',
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_blast_mixed(self) -> None:
        """
        Tests sequence typing using BLAST with a mixed scheme (DNA & peptide loci).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta),
            '--scheme-dir', str(self.input_db_mixed),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'blast',
            '--threads', '8',
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_srst2_nucl(self) -> None:
        """
        Tests sequence typing using SRST2 with a nucleotide scheme (including ST definitions).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(self.input_typing_reads[0]), str(self.input_typing_reads[1]),
            '--scheme-dir', str(self.input_db_nucl),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'srst2',
            '--trim-reads',
            '--srst2-max-unaligned-overlap', '123',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_srst2_mixed(self) -> None:
        """
        Tests sequence typing using SRST2 with a mixed scheme (DNA and peptide loci).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta),
            '--fastq-pe', str(self.input_typing_reads[0]), str(self.input_typing_reads[1]),
            '--scheme-dir', str(self.input_db_mixed),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'srst2',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
