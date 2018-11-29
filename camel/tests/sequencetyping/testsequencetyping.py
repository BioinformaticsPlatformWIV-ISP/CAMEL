import argparse
import unittest

import os
import tempfile

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.scripts.sequencetyping.mainsequencetyping import MainSequenceTyping


class TestSequenceTyping(unittest.TestCase):
    """
    Tests the sequence typing tool.
    """

    camel = Camel()
    running_dir = None

    # Input files
    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'])
    input_db_nucl = os.path.join(test_file_dir, 'typing', 'scheme_mlst_neisseria')
    input_db_protein = os.path.join(test_file_dir, 'typing', 'scheme_pora_neisseria')
    input_db_mixed = os.path.join(test_file_dir, 'typing', 'scheme_fhbp_neisseria')
    input_fasta = ToolIOFile(os.path.join(test_file_dir, 'typing', 'neisseria_mc58.fasta'))
    input_typing_reads = [
        ToolIOFile(os.path.join(test_file_dir, 'typing', 'S15BD05018_S58_L001_1.fastq')),
        ToolIOFile(os.path.join(test_file_dir, 'typing', 'S15BD05018_S58_L001_2.fastq'))
    ]

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(None, 'camel_', TestSequenceTyping.camel.config['temp_dir'])

    def test_typing_blast_nucl(self) -> None:
        """
        Tests sequence typing using BLAST with a nucleotide scheme (including ST definitions).
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            sample_name=None,
            fasta=self.input_fasta.path,
            fasta_name=os.path.basename(self.input_fasta.path),
            scheme_dir=self.input_db_nucl,
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            trim_reads=True,
            working_dir=self.running_dir,
            detection_method='blast',
            threads=8
        )
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_typing_blast_pept(self) -> None:
        """
        Tests sequence typing using BLAST with a peptide scheme.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            sample_name=None,
            fasta=self.input_fasta.path,
            fasta_name=os.path.basename(self.input_fasta.path),
            scheme_dir=self.input_db_protein,
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            trim_reads=True,
            working_dir=self.running_dir,
            detection_method='blast',
            threads=8
        )
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_typing_blast_mixed(self) -> None:
        """
        Tests sequence typing using BLAST with a mixed scheme (DNA & peptide loci).
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            sample_name=None,
            fasta=self.input_fasta.path,
            fasta_name=os.path.basename(self.input_fasta.path),
            scheme_dir=self.input_db_mixed,
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            trim_reads=True,
            working_dir=self.running_dir,
            detection_method='blast',
            threads=8
        )
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_typing_srst2_nucl(self) -> None:
        """
        Tests sequence typing using SRST2 with a nucleotide scheme (including ST definitions).
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            sample_name=None,
            fasta=None,
            fasta_name=None,
            fastq_pe=[x.path for x in self.input_typing_reads],
            fastq_pe_names=[os.path.basename(x.path) for x in self.input_typing_reads],
            scheme_dir=self.input_db_nucl,
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            trim_reads=True,
            working_dir=self.running_dir,
            detection_method='srst2',
            threads=8
        )
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_typing_srst2_mixed(self) -> None:
        """
        Tests sequence typing using SRST2 with a mixed scheme (DNA and peptide loci).
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            sample_name=None,
            fasta=self.input_fasta.path,
            fasta_name=os.path.basename(self.input_fasta.path),
            fastq_pe=[x.path for x in self.input_typing_reads],
            fastq_pe_names=[os.path.basename(x.path) for x in self.input_typing_reads],
            scheme_dir=self.input_db_mixed,
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            trim_reads=False,
            working_dir=self.running_dir,
            detection_method='srst2',
            threads=8
        )
        main = MainSequenceTyping(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)


if __name__ == '__main__':
    unittest.main()
