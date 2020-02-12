import unittest
from pathlib import Path

import tempfile

from camel.app.camel import Camel
from camel.scripts.genedetection.maingenedetection import MainGeneDetection


class TestGeneDetection(unittest.TestCase):
    """
    Tests the gene detection workflow.
    """

    camel = Camel()
    running_dir = None

    # Input files
    test_file_dir = Path(camel.config['testing']['testfiles_dir'])
    input_fasta = test_file_dir / 'workflows' / 'NC_002695.1.fasta'
    input_fasta_galaxy = test_file_dir / 'workflows' / 'dataset_12.dat'
    input_reads_no_hit = [test_file_dir / 'workflows' / 'ecoli_1.fastq',
                          test_file_dir / 'workflows' / 'ecoli_2.fastq']
    input_reads_raw = [test_file_dir / 'gene_detection' / 'reads-ds_1P.fastq',
                       test_file_dir / 'gene_detection' / 'reads-ds_2P.fastq']
    input_reads_raw_galaxy = [test_file_dir / 'workflows' / 'dataset_fwd_11.dat',
                              test_file_dir / 'workflows' / 'dataset_rev_10.dat']
    input_gene_detection_db = test_file_dir / 'gene_detection' / 'db'

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = Path(tempfile.mkdtemp(None, 'camel_', TestGeneDetection.camel.config['temp_dir']))

    def test_gene_detection_blast_fasta_input(self) -> None:
        """
        Tests the gene detection main script using blast.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(TestGeneDetection.input_fasta),
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir)
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_fasta_input_galaxy(self) -> None:
        """
        Tests the gene detection main script using blast with an input file in the Galaxy name style..
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(TestGeneDetection.input_fasta_galaxy),
            '--fasta-name', TestGeneDetection.input_fasta.name,
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir)
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_fasta_input_spaces(self) -> None:
        """
        Tests the gene detection main script using blast.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fasta', str(TestGeneDetection.input_fasta_galaxy),
            '--fasta-name', '"my reference genome.fasta"',
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir)
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_fastq_input(self) -> None:
        """
        Tests the gene detection main script using blast.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_raw_galaxy],
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--trim-reads'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_fastq_input_no_trim(self) -> None:
        """
        Tests the gene detection main script using blast.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_raw_galaxy],
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir)
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_galaxy_trimmomatic(self) -> None:
        """
        Tests the gene detection main script using blast with output generated by trimmomatic in Galaxy.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_raw_galaxy],
            '--fastq-pe-names',
            'Trimmomatic on Neisseria_2.fastq (R1 paired)', 'Trimmomatic on Neisseria_2.fastq (R2 paired)',
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir)
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_blast_invalid_filenames(self) -> None:
        """
        Tests the gene detection workflow with invalid filenames.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_raw_galaxy],
            '--fastq-pe-names', 'InvalidFilename 123', 'InvalidReverseReads 555',
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir)
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_srst2(self) -> None:
        """
        Tests the gene detection main script using SRST2.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_raw],
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'srst2',
            '--trim-reads'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_srst2_galaxy(self) -> None:
        """
        Tests the gene detection main script using SRST2 with Galaxy style inputs.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_raw_galaxy],
            '--fastq-pe-names', *[x.name for x in TestGeneDetection.input_reads_raw],
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'srst2'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_srst2_no_trim(self) -> None:
        """
        Tests the gene detection main script using SRST2 with Galaxy style inputs.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_raw],
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'srst2',
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_srst2_no_hits(self) -> None:
        """
        Tests the gene detection main script using SRST2 where no hits are detected.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_no_hit],
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'srst2',
            '--trim-reads'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_srst2_galaxy_trimmomatic(self) -> None:
        """
        Tests the gene detection main script using blast with output generated by trimmomatic in Galaxy.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_raw_galaxy],
            '--fastq-pe-names',
            'Trimmomatic on Neisseria_2.fastq (R1 paired)', 'Trimmomatic on Neisseria_2.fastq (R2 paired)',
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'srst2',
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)

    def test_gene_detection_kma(self) -> None:
        """
        Tests the gene detection main script using KMA.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        args = [
            '--fastq-pe', *[str(x) for x in TestGeneDetection.input_reads_raw],
            '--database-dir', str(TestGeneDetection.input_gene_detection_db),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--detection-method', 'kma',
            '--trim-reads'
        ]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
