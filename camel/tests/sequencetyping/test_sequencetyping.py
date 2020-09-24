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
    input_typing_reads = {
        'illumina': [test_file_dir / 'S15BD05018_S58_L001_1.fastq', test_file_dir / 'S15BD05018_S58_L001_2.fastq'],
        'iontorrent': ['/temp/bebog/dataset_typing/reads_neisseria.fastq'],
        'nanopore': [Path(
            '/data/testdata/neisseria/ERR2259087.fastq.gz')]
    }

    def test_typing_illumina_blast_nucl(self) -> None:
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

    def test_typing_illumina_blast_peptide(self) -> None:
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

    def test_typing_illumina_blast_mixed(self) -> None:
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

    def test_typing_illumina_srst2_nucl(self) -> None:
        """
        Tests sequence typing using SRST2 with a nucleotide scheme (including ST definitions).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(self.input_typing_reads['illumina'][0]), str(self.input_typing_reads['illumina'][1]),
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

    def test_typing_illumina_srst2_mixed(self) -> None:
        """
        Tests sequence typing using SRST2 with a mixed scheme (DNA and peptide loci).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta),
            '--fastq-pe', str(self.input_typing_reads['illumina'][0]), str(self.input_typing_reads['illumina'][1]),
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

    def test_typing_illumina_kma_nucl(self) -> None:
        """
        Tests sequence typing using KMA with a nucleotide scheme (including ST definitions).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fastq-pe', str(self.input_typing_reads['illumina'][0]), str(self.input_typing_reads['illumina'][1]),
            '--scheme-dir', str(self.input_db_nucl),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_illumina_kma_mixed(self) -> None:
        """
        Tests sequence typing using KMA with a mixed scheme (including ST definitions).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta),
            '--fastq-pe', str(self.input_typing_reads['illumina'][0]), str(self.input_typing_reads['illumina'][1]),
            '--scheme-dir', str(self.input_db_mixed),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_iontorrent_kma_nucl(self) -> None:
        """
        Tests sequence typing using KMA with a nucleotide scheme (including ST definitions).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fastq-se', str(self.input_typing_reads['iontorrent'][0]),
            '--scheme-dir', str(self.input_db_nucl),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--read-type', 'iontorrent',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_nanopore_kma_nucl(self) -> None:
        """
        Tests sequence typing using KMA with a nucleotide scheme (including ST definitions).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fastq-se', str(self.input_typing_reads['nanopore'][0]),
            '--scheme-dir', str(self.input_db_nucl),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--read-type', 'nanopore',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_iontorrent_kma_nucl_trim(self) -> None:
        """
        Tests sequence typing using KMA with a nucleotide scheme (including ST definitions).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fastq-se', str(self.input_typing_reads['iontorrent'][0]),
            '--scheme-dir', str(self.input_db_nucl),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--read-type', 'iontorrent',
            '--threads', '8',
            '--trim-reads'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_iontorrent_srst2_nucl(self) -> None:
        """
        Tests sequence typing using KMA with a nucleotide scheme (including ST definitions).
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fastq-se', str(self.input_typing_reads['iontorrent'][0]),
            '--scheme-dir', str(self.input_db_nucl),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'srst2',
            '--read-type', 'iontorrent',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_iontorrent_kma_mixed(self) -> None:
        """
        Tests sequence typing using KMA with a mixed scheme.
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta),
            '--fastq-se', str(self.input_typing_reads['iontorrent'][0]),
            '--scheme-dir', str(self.input_db_mixed),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--read-type', 'iontorrent',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)

    def test_typing_iontorrent_srst2_mixed(self) -> None:
        """
        Tests sequence typing using KMA with a mixed scheme.
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = [
            '--fasta', str(self.input_fasta),
            '--fastq-se', str(self.input_typing_reads['iontorrent'][0]),
            '--scheme-dir', str(self.input_db_mixed),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'srst2',
            '--read-type', 'iontorrent',
            '--threads', '8'
        ]
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(output_file_report.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
