import argparse
import unittest

import os
import tempfile

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.scripts.genedetection.maingenedetection import MainGeneDetection


class TestGeneDetection(unittest.TestCase):
    """
    Tests the gene detection workflow.
    """

    camel = Camel()
    running_dir = None

    # Input files
    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'])
    input_fasta = ToolIOFile(os.path.join(test_file_dir, 'workflows', 'NC_002695.1.fasta'))
    input_fasta_galaxy = ToolIOFile(os.path.join(test_file_dir, 'workflows', 'dataset_12.dat'))
    input_reads_no_hit = [ToolIOFile(os.path.join(test_file_dir, 'workflows', 'ecoli_1.fastq')),
                          ToolIOFile(os.path.join(test_file_dir, 'workflows', 'ecoli_2.fastq'))]
    input_reads_raw = [ToolIOFile('/data/camel/testdata/gene_detection/reads-ds_1P.fastq'),
                       ToolIOFile('/data/camel/testdata/gene_detection/reads-ds_2P.fastq')]
    input_reads_raw_galaxy = [ToolIOFile('/data/camel/testdata/workflows/dataset_fwd_11.dat'),
                              ToolIOFile('/data/camel/testdata/workflows/dataset_rev_10.dat')]
    input_gene_detection_db = os.path.join(test_file_dir, 'gene_detection', 'db')

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(None, 'camel_', TestGeneDetection.camel.config['temp_dir'])

    def __get_basic_arguments(self, report_path: str, detection_method: str) -> argparse.Namespace:
        """
        Returns the basic arguments for the main script.
        :param report_path: Report path
        :param detection_method: Detection method
        :return: Arguments
        """
        return argparse.Namespace(
            sample_name=None,
            fasta=None,
            fasta_name=None,
            fastq_pe=None,
            fastq_pe_names=None,
            database_dir=TestGeneDetection.input_gene_detection_db,
            output_html=report_path,
            output_dir=os.path.dirname(report_path),
            trim_reads=True,
            blast_min_percent_identity=90,
            blast_min_percent_coverage=70,
            srst2_min_cov=60,
            srst2_max_div=50,
            srst2_max_unaligned_overlap=150,
            srst2_max_mismatch=5,
            working_dir=self.running_dir,
            detection_method=detection_method,
            threads=4,
            report_include_fastq=False,
            kmers=55
        )

    def test_gene_detection_blast_fasta_input(self) -> None:
        """
        Tests the gene detection main script using blast.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = self.__get_basic_arguments(output_file_report, 'blast')
        args.fasta = TestGeneDetection.input_fasta.path
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_gene_detection_blast_fasta_input_galaxy(self) -> None:
        """
        Tests the gene detection main script using blast with an input file in the Galaxy name style..
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = self.__get_basic_arguments(output_file_report, 'blast')
        args.fasta = TestGeneDetection.input_fasta_galaxy.path
        args.fasta_name = TestGeneDetection.input_fasta.basename
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_gene_detection_blast_fasta_input_spaces(self) -> None:
        """
        Tests the gene detection main script using blast.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = self.__get_basic_arguments(output_file_report, 'blast')
        args.fasta = TestGeneDetection.input_fasta.path
        args.fasta_name = 'my reference genome.fasta'
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_gene_detection_blast_fastq_input(self) -> None:
        """
        Tests the gene detection main script using blast.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = self.__get_basic_arguments(output_file_report, 'blast')
        args.fastq_pe = [f.path for f in TestGeneDetection.input_reads_raw_galaxy]
        args.fastq_pe_names = [f.basename for f in TestGeneDetection.input_reads_raw]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_gene_detection_blast_fastq_input_no_trim(self) -> None:
        """
        Tests the gene detection main script using blast.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = self.__get_basic_arguments(output_file_report, 'blast')
        args.fastq_pe = [f.path for f in TestGeneDetection.input_reads_raw_galaxy]
        args.fastq_pe_names = [f.basename for f in TestGeneDetection.input_reads_raw]
        args.trim_reads = False
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_gene_detection_blast_galaxy_trimmomatic(self) -> None:
        """
        Tests the gene detection main script using blast with output generated by trimmomatic in Galaxy.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = self.__get_basic_arguments(output_file_report, 'blast')
        args.fastq_pe = [f.path for f in TestGeneDetection.input_reads_raw_galaxy]
        args.fastq_pe_names = ['Trimmomatic on Neisseria_2.fastq (R1 paired)',
                               'Trimmomatic on Neisseria_2.fastq (R2 paired)']
        args.trim_reads = False
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_gene_detection_blast_invalid_filenames(self) -> None:
        """
        Tests the gene detection workflow with invalid filenames.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = self.__get_basic_arguments(output_file_report, 'blast')
        args.fastq_pe = [f.path for f in TestGeneDetection.input_reads_raw_galaxy]
        args.fastq_pe_names = ['InvalidFilename 123', 'InvalidReverseReads 555']
        args.trim_reads = False
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_gene_detection_srst2(self) -> None:
        """
        Tests the gene detection main script using SRST2.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = self.__get_basic_arguments(output_file_report, 'srst2')
        args.fastq_pe = [f.path for f in TestGeneDetection.input_reads_raw]
        args.fastq_pe_names = [os.path.basename(p.path) for p in TestGeneDetection.input_reads_raw]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_gene_detection_srst2_galaxy(self) -> None:
        """
        Tests the gene detection main script using SRST2 with Galaxy style inputs.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = self.__get_basic_arguments(output_file_report, 'srst2')
        args.fastq_pe = [f.path for f in TestGeneDetection.input_reads_raw_galaxy]
        args.fastq_pe_names = [os.path.basename(p.path) for p in TestGeneDetection.input_reads_raw]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_gene_detection_srst2_no_trim(self) -> None:
        """
        Tests the gene detection main script using SRST2 with Galaxy style inputs.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = self.__get_basic_arguments(output_file_report, 'srst2')
        args.trim_reads = None
        args.fastq_pe = [f.path for f in TestGeneDetection.input_reads_raw_galaxy]
        args.fastq_pe_names = [os.path.basename(p.path) for p in TestGeneDetection.input_reads_raw]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_gene_detection_srst2_no_hits(self) -> None:
        """
        Tests the gene detection main script using SRST2 where no hits are detected.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = self.__get_basic_arguments(output_file_report, 'srst2')
        args.fastq_pe = [f.path for f in TestGeneDetection.input_reads_no_hit]
        args.fastq_pe_names = [os.path.basename(p.path) for p in TestGeneDetection.input_reads_no_hit]
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_gene_detection_srst2_galaxy_trimmomatic(self) -> None:
        """
        Tests the gene detection main script using blast with output generated by trimmomatic in Galaxy.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = self.__get_basic_arguments(output_file_report, 'srst2')
        args.fastq_pe = [f.path for f in TestGeneDetection.input_reads_raw_galaxy]
        args.fastq_pe_names = ['Trimmomatic on Neisseria_2.fastq (R1 paired)',
                               'Trimmomatic on Neisseria_2.fastq (R2 paired)']
        args.trim_reads = False
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)


if __name__ == '__main__':
    unittest.main()
