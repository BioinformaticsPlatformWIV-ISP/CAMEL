import argparse
import unittest

from camel.app.camel import Camel
from camel.app.components.phylogeny.snpphylogenyutils import SnpPhylogenyUtils
from camel.app.error.invalidinputerror import InvalidInputError


class TestSnpPhylogenyUtils(unittest.TestCase):
    """
    Tests the SRST2 base tools.
    """
    camel = Camel()

    def test_valid_samples(self):
        """
        Tests if valid samples are correctly added.
        :return: None
        """
        args = argparse.Namespace(
            sample=[
                ('S1', 'S1_1.fastq', '/path/to/rS1_1.fastq', 'S1_2.fastq', '/path/to/S1_2.fastq'),
                ('S2', 'S2_1.fastq', '/path/to/rS2_1.fastq', 'S2_2.fastq', '/path/to/S2_2.fastq'),
                ('S3', 'S3_1.fastq', '/path/to/rS3_1.fastq', 'S3_2.fastq', '/path/to/S3_2.fastq'),
                ('S4', 'S4_1.fastq', '/path/to/rS4_1.fastq', 'S4_2.fastq', '/path/to/S4_2.fastq')
            ])
        samples = SnpPhylogenyUtils.extract_samples(args)
        self.assertEqual(len(samples), 4)

    def test_valid_samples_from_reads(self):
        """
        Tests if valid samples are correctly added from the read names.
        :return: None
        """
        args = argparse.Namespace(
            sample=[
                ('', 'S1_1.fastq', '/path/to/rS1_1.fastq', 'S1_2.fastq', '/path/to/S1_2.fastq'),
                ('', 'S2_1.fastq', '/path/to/rS2_1.fastq', 'S2_2.fastq', '/path/to/S2_2.fastq'),
                ('', 'S3_1.fastq', '/path/to/rS3_1.fastq', 'S3_2.fastq', '/path/to/S3_2.fastq'),
                ('', 'S4_1.fastq', '/path/to/rS4_1.fastq', 'S4_2.fastq', '/path/to/S4_2.fastq')
            ])
        samples = SnpPhylogenyUtils.extract_samples(args)
        self.assertEqual(len(samples), 4)

    def test_duplicate_samples(self):
        """
        Tests if duplicate samples raise an error.
        :return: None
        """
        args = argparse.Namespace(
            sample=[
                ('S1', 'S1_1.fastq', '/path/to/rS1_1.fastq', 'S1_2.fastq', '/path/to/S1_2.fastq'),
                ('S2', 'S2_1.fastq', '/path/to/rS2_1.fastq', 'S2_2.fastq', '/path/to/S2_2.fastq'),
                ('S3', 'S3_1.fastq', '/path/to/rS3_1.fastq', 'S3_2.fastq', '/path/to/S3_2.fastq'),
                ('S3', 'S3_1.fastq', '/path/to/rS3_1.fastq', 'S3_2.fastq', '/path/to/S3_2.fastq'),
                ('S4', 'S4_1.fastq', '/path/to/rS4_1.fastq', 'S4_2.fastq', '/path/to/S4_2.fastq')
            ])
        with self.assertRaises(InvalidInputError):
            SnpPhylogenyUtils.extract_samples(args)
