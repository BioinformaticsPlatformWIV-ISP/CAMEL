import argparse
import unittest
from typing import List

import os
import tempfile

from camel.app.camel import Camel
from camel.scripts.snpphylogeny.maincfsanphylo import MainCfsanPhylo
from camel.scripts.snpphylogeny.mainsamtoolsphylo import MainSamtoolsPhylo


class TestSnpPhylogenyPipelines(unittest.TestCase):
    """
    Tests the SNP phylogeny pipelines.
    """
    camel = Camel.get_instance()
    running_dir = None
    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'])
    reference_fasta = os.path.join(
        camel.config['testing']['testfiles_dir'], 'phylogeny', 'lm_subset', 'reference.fasta')

    def setUp(self) -> None:
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(prefix='camel_', dir=TestSnpPhylogenyPipelines.camel.config['temp_dir'])

    @staticmethod
    def __get_samples(gzipped: bool = False) -> List[List[str]]:
        """
        Returns the input samples.
        :return: Input samples
        """
        reads_dir = os.path.join(TestSnpPhylogenyPipelines.test_file_dir, 'phylogeny', 'lm_subset')
        return [[
                f'Sample_{i}',
                f"lm{i}_1.fastq{'.gz' if gzipped else ''}",
                os.path.join(reads_dir, f"lm{i}_1.fastq{'.gz' if gzipped else ''}"),
                f"lm{i}_2.fastq{'.gz' if gzipped else ''}",
                os.path.join(reads_dir, f"lm{i}_2.fastq{'.gz' if gzipped else ''}")] for i in range(1, 5)]

    def test_samtools_phylogeny(self) -> None:
        """
        Tests the Samtools Phylogeny pipeline.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            sample=TestSnpPhylogenyPipelines.__get_samples(),
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            working_dir=self.running_dir,
            reference_name=os.path.basename(TestSnpPhylogenyPipelines.reference_fasta),
            reference=TestSnpPhylogenyPipelines.reference_fasta,
            trim_reads=True,
            ploidy='diploid',
            calling_method='multiallelic',
            min_total_depth=5,
            min_forward_depth=0,
            min_reverse_depth=0,
            min_snp_quality=20,
            min_mapping_quality=25,
            min_distance=3,
            keep_best=False,
            min_zscore=0.5,
            y_mult=3,
            missing_data='complete_deletion',
            branch_swap='none',
            site_cov_cutoff=50,
            bootstraps=10,
            ml_method='spr3',
            threads=8,
            report_include_bam=False,
            include_ref=True
        )
        main = MainSamtoolsPhylo(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_cfsan_phylogeny(self) -> None:
        """
        Tests the CFSAN Phylogeny pipeline.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            sample=TestSnpPhylogenyPipelines.__get_samples(),
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            working_dir=self.running_dir,
            reference_name=os.path.basename(TestSnpPhylogenyPipelines.reference_fasta),
            reference=TestSnpPhylogenyPipelines.reference_fasta,
            trim_reads=True,
            missing_data='complete_deletion',
            branch_swap='none',
            site_cov_cutoff=50,
            bootstraps=10,
            ml_method='spr3',
            threads=8,
            export_bam=False,
            selected_matrix='regular',
            report_include_bam=False
        )
        main = MainCfsanPhylo(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_cfsan_phylogeny_gz_no_trim(self) -> None:
        """
        Tests the CFSAN Phylogeny pipeline with gzipped input files and no trimming enabled.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            sample=TestSnpPhylogenyPipelines.__get_samples(True),
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            working_dir=self.running_dir,
            reference_name=os.path.basename(TestSnpPhylogenyPipelines.reference_fasta),
            reference=TestSnpPhylogenyPipelines.reference_fasta,
            trim_reads=False,
            missing_data='complete_deletion',
            branch_swap='none',
            site_cov_cutoff=50,
            bootstraps=10,
            ml_method='spr3',
            threads=8,
            export_bam=False,
            selected_matrix='regular',
            report_include_bam=False
        )
        main = MainCfsanPhylo(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)
