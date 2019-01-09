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
    input_reads_no_hit = [ToolIOFile(os.path.join(test_file_dir, 'workflows', 'ecoli_1.fastq')),
                          ToolIOFile(os.path.join(test_file_dir, 'workflows', 'ecoli_2.fastq'))]
    input_reads_raw = [ToolIOFile('/data/camel/testdata/gene_detection/reads-ds_1P.fastq'),
                       ToolIOFile('/data/camel/testdata/gene_detection/reads-ds_2P.fastq')]
    input_gene_detection_db = os.path.join(test_file_dir, 'gene_detection', 'db')

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(None, 'camel_', TestGeneDetection.camel.config['temp_dir'])

    def test_gene_detection_blast(self) -> None:
        """
        Tests the gene detection main script using blast.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            sample_name=None,
            fasta=self.input_fasta.path,
            fasta_name=os.path.basename(self.input_fasta.path),
            fastq_pe=None,
            fastq_pe_names=None,
            database_dir=self.input_gene_detection_db,
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            trim_reads=True,
            blast_min_percent_identity=90,
            blast_min_percent_coverage=70,
            working_dir=self.running_dir,
            detection_method='blast',
            threads=4
        )
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_gene_detection_srst2(self) -> None:
        """
        Tests the gene detection main script using SRST2.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            sample_name=None,
            fasta_name=None,
            fasta=None,
            fastq_pe=[f.path for f in TestGeneDetection.input_reads_raw],
            fastq_pe_names=[os.path.basename(p.path) for p in TestGeneDetection.input_reads_raw],
            database_dir=self.input_gene_detection_db,
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            trim_reads=True,
            srst2_min_cov=50,
            srst2_max_div=30,
            working_dir=self.running_dir,
            detection_method='srst2',
            threads=4,
            report_include_fastq=False
        )
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_gene_detection_srst2_no_hits(self) -> None:
        """
        Tests the gene detection main script using SRST2 where no hits are detected.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            sample_name=None,
            fasta_name=None,
            fasta=None,
            fastq_pe=[f.path for f in TestGeneDetection.input_reads_no_hit],
            fastq_pe_names=[os.path.basename(p.path) for p in TestGeneDetection.input_reads_no_hit],
            database_dir=self.input_gene_detection_db,
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            trim_reads=True,
            srst2_min_cov=50,
            srst2_max_div=30,
            working_dir=self.running_dir,
            detection_method='srst2',
            threads=4,
            report_include_fastq=False
        )
        main = MainGeneDetection(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)


if __name__ == '__main__':
    unittest.main()
