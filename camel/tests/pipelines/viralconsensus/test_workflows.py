import logging
import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.components.vcf.vcfutils import VCFUtils
from camel.app.components.workflows.trimmingilluminawrapper import TrimmingIlluminaWrapper
from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.io.tooliofile import ToolIOFile
from camel.scripts.viralconsensuspipeline.workflows.applyvariants import ApplyVariants
from camel.scripts.viralconsensuspipeline.workflows.callvariants import CallVariants
from camel.scripts.viralconsensuspipeline.workflows.filtervariants import FilterVariants
from camel.scripts.viralconsensuspipeline.workflows.readmappingworkflow import ReadMappingWorkflow


class TestWorkflows(CamelTestSuite):
    """
    Tests for the sub-workflows for the viral consensus pipeline.
    """

    def test_read_mapping_workflow_illumina(self) -> None:
        """
        Tests the read mapping workflow.
        :return: None
        """
        test_file_dir_illumina = CamelTestSuite.get_test_file_dir('read_mapping')
        fasta_ref = test_file_dir_illumina / 'reference.fasta'
        fq_fwd = test_file_dir_illumina / 'reads_1.fastq'
        fq_rev = test_file_dir_illumina / 'reads_2.fastq'
        fastq_input = FastqInput('illumina', [ToolIOFile(fq_fwd), ToolIOFile(fq_rev)], is_pe=True)
        workflow = ReadMappingWorkflow(self.running_dir)
        output = workflow.run(fastq_input, fasta_ref)
        self.assertTrue(output.path_bam.exists())
        self.assertGreater(output.path_bam.stat().st_size, 0)
        self.assertGreater(len(output.informs), 0)

    def test_read_mapping_workflow_illumina_depth_cutoff(self) -> None:
        """
        Tests the read mapping workflow.
        :return: None
        """
        test_file_dir_illumina = CamelTestSuite.get_test_file_dir('read_mapping')
        fasta_ref = test_file_dir_illumina / 'reference.fasta'
        fq_fwd = test_file_dir_illumina / 'reads_1.fastq'
        fq_rev = test_file_dir_illumina / 'reads_2.fastq'
        fastq_input = FastqInput('illumina', [ToolIOFile(fq_fwd), ToolIOFile(fq_rev)], is_pe=True)
        workflow = ReadMappingWorkflow(self.running_dir)
        output = workflow.run(fastq_input, fasta_ref, gap_depth_cutoff=1, gap_len_cutoff=3)
        self.assertTrue(output.path_bam.exists())
        self.assertGreater(output.path_bam.stat().st_size, 0)
        self.assertGreater(len(output.informs), 0)

    def test_read_mapping_workflow_illumina_trimmed(self) -> None:
        """
        Tests the read mapping workflow on trimmed reads.
        :return: None
        """
        # Input data
        test_file_dir_illumina = CamelTestSuite.get_test_file_dir('read_mapping')
        fasta_ref = test_file_dir_illumina / 'reference.fasta'
        fq_fwd = test_file_dir_illumina / 'reads_1.fastq'
        fq_rev = test_file_dir_illumina / 'reads_2.fastq'

        # Trim reads
        dir_trimming = self.running_dir / 'trim'
        wrapper = TrimmingIlluminaWrapper(dir_trimming)
        wrapper.run_workflow([fq_fwd, fq_rev], threads=4)
        fastq_input = FastqInput(
            read_type='illumina',
            pe=wrapper.output.trimmed_reads_pe,
            se_fwd=wrapper.output.trimmed_reads_se_fwd,
            se_rev=wrapper.output.trimmed_reads_se_rev,
            is_trimmed=True, is_pe=True)

        # Run read mapping workflow
        workflow = ReadMappingWorkflow(self.running_dir)
        output = workflow.run(fastq_input, fasta_ref, gap_len_cutoff=3, gap_depth_cutoff=12)
        self.assertTrue(output.path_bam.exists())
        self.assertGreater(output.path_bam.stat().st_size, 0)
        self.assertGreater(len(output.informs), 0)

    def test_read_mapping_workflow_ont(self) -> None:
        """
        Tests the read mapping workflow.
        :return: None
        """
        test_file_dir_ont = CamelTestSuite.get_test_file_dir('minion')
        fq_ont = test_file_dir_ont / 'fastq_minion_stec.fastq'
        fasta_ref_ont = test_file_dir_ont / 'NC_002695.2.fasta'
        fastq_input = FastqInput('nanopore', se=[ToolIOFile(fq_ont)], is_pe=False)
        workflow = ReadMappingWorkflow(self.running_dir)
        output = workflow.run(fastq_input, fasta_ref_ont)
        self.assertTrue(output.path_bam.exists())
        self.assertGreater(output.path_bam.stat().st_size, 0)
        self.assertGreater(len(output.informs), 0)

    def test_read_mapping_workflow_ont_no_depth(self) -> None:
        """
        Tests the read mapping workflow.
        :return: None
        """
        test_file_dir_ont = CamelTestSuite.get_test_file_dir('minion')
        fq_ont = test_file_dir_ont / 'fastq_minion_stec.fastq'
        fasta_ref_ont = test_file_dir_ont / 'H1N1_HA.fasta'
        fastq_input = FastqInput('nanopore', se=[ToolIOFile(fq_ont)], is_pe=False)
        workflow = ReadMappingWorkflow(self.running_dir)
        output = workflow.run(fastq_input, fasta_ref_ont)
        self.assertTrue(output.path_bam.exists())
        self.assertGreater(output.path_bam.stat().st_size, 0)
        self.assertGreater(len(output.informs), 0)

    def test_call_variants_bcftools_ont(self) -> None:
        """
        Runs the variant calling workflow with ONT data using bcftools.
        :return: None
        """
        dir_test = CamelTestSuite.get_test_file_dir('clair3')
        workflow = CallVariants(self.running_dir)
        output = workflow.run(dir_test / 'bsubtilis_ont.bam', dir_test / 'bsubtilis.fa', 'nanopore', 'bcftools', {})
        self.assertTrue(output.path_vcf.exists())
        self.assertGreater(VCFUtils.count_variants(output.path_vcf), 0)
        logging.info(f'Stats: {output.stats}')

    def test_call_variants_bcftools_illumina(self) -> None:
        """
        Runs the variant calling workflow with ONT data using bcftools.
        :return: None
        """
        dir_test = CamelTestSuite.get_test_file_dir('clair3')
        workflow = CallVariants(self.running_dir)
        output = workflow.run(
            dir_test / 'bsubtilis_illumina.bam', dir_test / 'bsubtilis.fa', 'illumina', 'bcftools', {})
        self.assertTrue(output.path_vcf.exists())
        self.assertGreater(VCFUtils.count_variants(output.path_vcf), 0)
        self.assertGreater(len(output.informs), 0)
        logging.info(f'Stats: {output.stats}')

    def test_call_variants_clair3_ont(self) -> None:
        """
        Runs the variant calling workflow with ONT data using bcftools.
        :return: None
        """
        dir_test = CamelTestSuite.get_test_file_dir('clair3')
        model = Path(self.camel.config['db_root']) / 'clair3' / 'models' / 'ont'
        workflow = CallVariants(self.running_dir)
        output = workflow.run(
            dir_test / 'bsubtilis_ont.bam', dir_test / 'bsubtilis.fa', 'nanopore', 'clair3', {'model': model},
            threads=8)
        self.assertTrue(output.path_vcf.exists())
        self.assertGreater(VCFUtils.count_variants(output.path_vcf), 0)
        self.assertGreater(len(output.informs), 0)
        logging.info(f'Stats: {output.stats}')

    def test_call_variants_clair3_illumina(self) -> None:
        """
        Runs the variant calling workflow with ONT data using bcftools.
        :return: None
        """
        dir_test = CamelTestSuite.get_test_file_dir('clair3')
        model = Path(self.camel.config['db_root']) / 'clair3' / 'models' / 'ilmn'
        workflow = CallVariants(self.running_dir)
        output = workflow.run(
            dir_test / 'bsubtilis_illumina.bam', dir_test / 'bsubtilis.fa', 'illumina', 'clair3', {'model': model},
            threads=4)
        self.assertTrue(output.path_vcf.exists())
        self.assertGreater(VCFUtils.count_variants(output.path_vcf), 0)
        self.assertGreater(len(output.informs), 0)
        logging.info(f'Stats: {output.stats}')

    def test_filter_variants_bcftools(self) -> None:
        """
        Tests the filtering of variants using bcftools.
        :return: None
        """
        dir_test = CamelTestSuite.get_test_file_dir('variant_calling')
        workflow = FilterVariants(self.running_dir)
        filters = {'min_af': 0.5, 'min_dp': 20, 'min_qual': 25}
        out = workflow.run(dir_test / 'unfiltered_variants-myco.vcf', 'bcftools', filters)
        self.assertTrue(out.path_vcf.exists())
        self.assertGreater(len(out.informs), 0)
        self.assertGreater(VCFUtils.count_variants(out.path_vcf), 0)

    def test_filter_variants_clair3(self) -> None:
        """
        Tests the filtering of variants using Clair3.
        :return: None
        """
        dir_test = CamelTestSuite.get_test_file_dir('variant_calling')
        workflow = FilterVariants(self.running_dir)
        filters = {'min_af': 0.5, 'min_dp': 20, 'min_qual': 25}
        out = workflow.run(dir_test / 'unfiltered_variants-myco-clair3.vcf.gz', 'clair3', filters)
        self.assertTrue(out.path_vcf.exists())
        self.assertGreater(len(out.informs), 0)
        self.assertGreater(VCFUtils.count_variants(out.path_vcf), 0)

    def test_apply_variants(self) -> None:
        """
        Tests the apply_variants workflow.
        :return: None
        """
        dir_test = CamelTestSuite.get_test_file_dir('variant_calling')
        workflow = ApplyVariants(self.running_dir)
        workflow.run(
            fasta_in=dir_test / 'ref-myco.fasta',
            vcf_in=dir_test / 'filtered_variants-myco.vcf.gz',
            name='updated'
        )


if __name__ == '__main__':
    unittest.main()
