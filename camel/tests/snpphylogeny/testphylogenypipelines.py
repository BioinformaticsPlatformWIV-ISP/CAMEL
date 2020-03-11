import unittest
from pathlib import Path
from typing import List

import os
import tempfile

from camel.app.camel import Camel
from camel.scripts.snpphylogeny.maincfsanphylo import MainCfsanPhylo
from camel.scripts.snpphylogeny.mainsamtoolsphylo import MainSamtoolsPhylo
from camel.tests import longRunningTest


class TestSnpPhylogenyPipelines(unittest.TestCase):
    """
    Tests the SNP phylogeny pipelines.
    """
    camel = Camel.get_instance()
    running_dir = None
    test_file_dir = Path(camel.config['testing']['testfiles_dir']) / 'phylogeny' / 'lm_subset'
    reference_fasta = test_file_dir / 'reference.fasta'

    def setUp(self) -> None:
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = Path(
            tempfile.mkdtemp(prefix='camel_', dir=TestSnpPhylogenyPipelines.camel.config['temp_dir']))

    @staticmethod
    def __get_samples(gzipped: bool = False) -> List[str]:
        """
        Returns the input samples.
        :return: Input samples
        """
        samples = []
        pattern = '*_1.fastq' + ('.gz' if gzipped else '')
        for index, fq_fwd in enumerate(TestSnpPhylogenyPipelines.test_file_dir.glob(pattern), 1):
            fq_rev = fq_fwd.parent / fq_fwd.name.replace('_1', '_2')
            samples.extend([
                '--sample', f'sample_{index}', fq_fwd.name, str(fq_fwd), fq_rev.name, str(fq_rev)])
        return samples

    @longRunningTest()
    def test_samtools_phylogeny(self) -> None:
        """
        Tests the Samtools Phylogeny pipeline.
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = TestSnpPhylogenyPipelines.__get_samples() + [
            '--reference', str(TestSnpPhylogenyPipelines.reference_fasta),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--trim-reads', '--include-ref',
            '--min-total-depth', '1',
            '--min-forward-depth', '0',
            '--min-reverse-depth', '0',
            '--min-distance', '2'
        ]
        main = MainSamtoolsPhylo(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    @longRunningTest()
    def test_cfsan_phylogeny(self) -> None:
        """
        Tests the CFSAN Phylogeny pipeline.
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = TestSnpPhylogenyPipelines.__get_samples() + [
            '--reference', str(TestSnpPhylogenyPipelines.reference_fasta),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--trim-reads',
            '--missing-data', 'complete_deletion',
            '--branch-swap', 'none',
            '--site-cov-cutoff', '50',
            '--bootstraps', '10'
        ]
        main = MainCfsanPhylo(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    @longRunningTest()
    def test_cfsan_phylogeny_gz_no_trim(self) -> None:
        """
        Tests the CFSAN Phylogeny pipeline with gzipped input files and no trimming enabled.
        :return: None
        """
        output_file_report = Path(self.running_dir) / 'report' / 'report.html'
        args = TestSnpPhylogenyPipelines.__get_samples(True) + [
            '--reference', str(TestSnpPhylogenyPipelines.reference_fasta),
            '--output-html', str(output_file_report),
            '--output-dir', str(output_file_report.parent),
            '--working-dir', str(self.running_dir),
            '--missing-data', 'complete_deletion',
            '--ml-method', 'spr3',
            '--site-cov-cutoff', '50',
            '--bootstraps', '10'
        ]
        main = MainCfsanPhylo(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)
