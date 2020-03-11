import argparse
import unittest

import os
import tempfile

from camel.app.camel import Camel
from camel.app.command.command import Command
from camel.app.tools.srst2.srst2gene import Srst2Gene
from camel.scripts.srst2.mainsrst2gene import MainSrst2Gene
from camel.scripts.srst2.mainsrst2mlst import MainSrst2Mlst
from camel.tests import longRunningTest


class TestSRST2(unittest.TestCase):
    """
    Tests the SRST2 base tools.
    """
    camel = Camel()
    running_dir = None

    # Input files
    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'])
    input_pe_reads = [os.path.join(test_file_dir, 'srst2', f'ERR178148-ds_1.fastq'),
                      os.path.join(test_file_dir, 'srst2', f'ERR178148-ds_2.fastq')]
    input_gene_db_path = os.path.join(test_file_dir, 'srst2', 'resfinder_clustered_srst2.fasta')
    input_mlst_db_fasta = os.path.join(test_file_dir, 'srst2', 'mlst_warwick_all_seq.fasta')
    input_mlst_db_profiles = os.path.join(test_file_dir, 'srst2', 'profiles.tsv')

    def setUp(self) -> None:
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(prefix='camel_', dir=TestSRST2.camel.config['temp_dir'])

    def test_dependencies(self) -> None:
        """
        Tests if the tool dependencies are available.
        :return: None
        """
        srst2 = Srst2Gene(Camel.get_instance())
        for dependency in srst2.dependencies:
            command = Command(f'module load {dependency};')
            command.run_command(self.running_dir)
            self.assertEqual(command.returncode, 0, f"Dependency '{dependency}' cannot be loaded")

    def test_pe_gene_detection(self) -> None:
        """
        Tests gene detection using PE reads as input
        :return: None
        """
        output_file = os.path.join(self.running_dir, 'srst2_tab_out.tsv')
        args = argparse.Namespace(
            fastq_pe=TestSRST2.input_pe_reads,
            gene_fasta=TestSRST2.input_gene_db_path,
            gene_fasta_name=os.path.basename(TestSRST2.input_gene_db_path),
            working_dir=self.running_dir,
            output_tsv=output_file,
            min_coverage=50,
            max_divergence=20,
            min_depth=None,
            min_edge_depth=None,
            max_unaligned_overlap=None)
        srst2_gene = MainSrst2Gene(args)
        srst2_gene.run()
        self.assertGreater(os.path.getsize(output_file), 0)

    def test_se_gene_detection(self) -> None:
        """
        Tests gene detection using SE reads as input
        :return: None
        """
        output_file = os.path.join(self.running_dir, 'srst2_tab_out.tsv')
        args = argparse.Namespace(
            fastq_se=TestSRST2.input_pe_reads[0],
            fastq_pe=None,
            gene_fasta=TestSRST2.input_gene_db_path,
            gene_fasta_name=os.path.basename(TestSRST2.input_gene_db_path),
            working_dir=self.running_dir,
            output_tsv=output_file,
            min_coverage=50,
            max_divergence=20,
            min_depth=None,
            min_edge_depth=None,
            max_unaligned_overlap=None)
        srst2_gene = MainSrst2Gene(args)
        srst2_gene.run()
        self.assertGreater(os.path.getsize(output_file), 0)

    @longRunningTest()
    def test_se_typing(self) -> None:
        """
        Tests typing using SE input.
        :return: None
        """
        output_file = os.path.join(self.running_dir, 'srst2_tab_out.tsv')
        args = argparse.Namespace(
            fastq_se=TestSRST2.input_pe_reads[0],
            fastq_pe=None,
            locus_fasta=TestSRST2.input_mlst_db_fasta,
            profiles_tsv=TestSRST2.input_mlst_db_profiles,
            working_dir=self.running_dir,
            output_tsv=output_file,
            max_mismatch=4,
            min_depth=None,
            min_edge_depth=None,
            max_unaligned_overlap=None)
        srst2_mlst = MainSrst2Mlst(args)
        srst2_mlst.run()
        self.assertGreater(os.path.getsize(output_file), 0)

    @longRunningTest()
    def test_pe_typing(self) -> None:
        """
        Tests typing using PE input without profiles
        :return: None
        """
        output_file = os.path.join(self.running_dir, 'srst2_tab_out.tsv')
        args = argparse.Namespace(
            fastq_pe=TestSRST2.input_pe_reads,
            locus_fasta=TestSRST2.input_mlst_db_fasta,
            profiles_tsv=None,
            working_dir=self.running_dir,
            output_tsv=output_file,
            max_mismatch=4,
            min_depth=None,
            min_edge_depth=None,
            max_unaligned_overlap=None)
        srst2_mlst = MainSrst2Mlst(args)
        srst2_mlst.run()
        self.assertGreater(os.path.getsize(output_file), 0)
