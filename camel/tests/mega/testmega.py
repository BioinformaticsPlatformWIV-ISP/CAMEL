import argparse
import unittest

import os
import tempfile

from camel.app.camel import Camel
from camel.scripts.mega.mainmega import MainMega


class TestMEGA(unittest.TestCase):
    """
    Tests the MEGA tool.
    """
    camel = Camel()
    running_dir = None

    # Input files
    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'])
    input_snp_matrix = os.path.join(test_file_dir, 'mega', 'sample_snp_matrix.fasta')
    input_vcf_files = [
        os.path.join(test_file_dir, 'mega', 'variants-s1.vcf'),
        os.path.join(test_file_dir, 'mega', 'variants-s2.vcf'),
        os.path.join(test_file_dir, 'mega', 'variants-s3.vcf'),
        os.path.join(test_file_dir, 'mega', 'variants-s4.vcf')
    ]

    def setUp(self) -> None:
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(prefix='camel_', dir=TestMEGA.camel.config['temp_dir'])

    def test_model_selection(self) -> None:
        """
        Tests the model selection
        :return: None
        """
        output_file = os.path.join(self.running_dir, 'output_model.tsv')
        args = argparse.Namespace(
            fasta=TestMEGA.input_snp_matrix,
            working_dir=self.running_dir,
            action='model',
            output_model=output_file,
            missing_data='complete_deletion',
            site_cov_cutoff=50,
            branch_swap='weak',
            threads=4
        )
        mega = MainMega(args)
        mega.run()
        self.assertGreater(os.path.getsize(output_file), 0)

    def test_tree_building(self) -> None:
        """
        Tests the tree building.
        :return: None
        """
        output_file = os.path.join(self.running_dir, 'output_tree.nwk')
        args = argparse.Namespace(
            fasta=TestMEGA.input_snp_matrix,
            working_dir=self.running_dir,
            action='tree',
            model='T92',
            rates='G',
            output_tree=output_file,
            missing_data='use_all_sites',
            site_cov_cutoff=50,
            branch_swap='moderate',
            ml_method='spr3',
            bootstraps=10,
            threads=4
        )
        mega = MainMega(args)
        mega.run()
        self.assertGreater(os.path.getsize(output_file), 0)

    def test_both(self) -> None:
        """
        Tests model selection + tree building.
        :return: None
        """
        output_file_tree = os.path.join(self.running_dir, 'output_tree.nwk')
        output_file_model = os.path.join(self.running_dir, 'output_model.tsv')
        args = argparse.Namespace(
            fasta=TestMEGA.input_snp_matrix,
            working_dir=self.running_dir,
            action='both',
            output_tree=output_file_tree,
            output_model=output_file_model,
            missing_data='partial_deletion',
            site_cov_cutoff=50,
            branch_swap='weak',
            ml_method='nni',
            bootstraps=10,
            threads=4
        )
        mega = MainMega(args)
        mega.run()
        self.assertGreater(os.path.getsize(output_file_model), 0)
        self.assertGreater(os.path.getsize(output_file_tree), 0)

    def test_export_snp_matrix(self) -> None:
        """
        Tests the export SNP matrix function (starting from multiple VCF input files).
        :return: None
        """
        output_file_snp_matrix = os.path.join(self.running_dir, 'snps.fasta')
        output_file_model = os.path.join(self.running_dir, 'output_model.tsv')
        args = argparse.Namespace(
            vcf=[[f, os.path.basename(f)] for f in TestMEGA.input_vcf_files],
            fasta=None,
            working_dir=self.running_dir,
            action='model',
            output_model=output_file_model,
            output_snp_matrix=output_file_snp_matrix,
            missing_data='use_all_sites',
            site_cov_cutoff=50,
            branch_swap='very_strong',
            threads=4
        )
        mega = MainMega(args)
        mega.run()
        self.assertGreater(os.path.getsize(output_file_model), 0)
        self.assertGreater(os.path.getsize(output_file_snp_matrix), 0)
