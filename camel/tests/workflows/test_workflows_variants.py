import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.components.workflows.variantcallingwrapper import VariantCallingWrapper
from camel.app.components.workflows.variantfilteringwrapper import VariantFilteringWrapper
from camel.app.io.tooliofile import ToolIOFile


class TestWorkflowsVariants(CamelTestSuite):
    """
    Tests the Snakemake variant calling related workflows.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir()

    # Variant calling
    input_fastq_pe = [
        test_file_dir / 'variant_calling' / 'toy' / 'r_1.fq',
        test_file_dir / 'variant_calling' / 'toy' / 'r_2.fq']
    input_fasta_ref = test_file_dir / 'variant_calling' / 'toy' / 'toy.fasta'

    # Variant filtering
    input_vcf = test_file_dir / 'variant_calling' / 'unfiltered_variants-myco.vcf.gz'
    input_bam = test_file_dir / 'variant_calling' / 'alignment.bam'

    def test_variant_calling_illumina(self) -> None:
        """
        Tests the variant calling workflow.
        :return: None
        """
        wrapper = VariantCallingWrapper(self.running_dir)
        fastq_input = FastqInput('illumina', pe=[ToolIOFile(x) for x in TestWorkflowsVariants.input_fastq_pe])
        ref_info = {
            'name': TestWorkflowsVariants.input_fasta_ref.name,
            'path': str(TestWorkflowsVariants.input_fasta_ref)
        }
        wrapper.run_workflow(ref_info, 'test_sample', fastq_input, {}, 1)
        self.assertGreater(wrapper.output.bam_file.size, 0)
        self.assertGreater(wrapper.output.vcf_unfiltered.size, 0)

    def test_variant_calling_single_end(self) -> None:
        """
        Tests the variant calling workflow with IonTorrent input data.
        :return: None
        """
        wrapper = VariantCallingWrapper(self.running_dir)
        fastq_input = FastqInput('iontorrent', se=[ToolIOFile(TestWorkflowsVariants.input_fastq_pe[0])], is_pe=False)
        ref_info = {
            'name': TestWorkflowsVariants.input_fasta_ref.name,
            'path': str(TestWorkflowsVariants.input_fasta_ref)
        }
        wrapper.run_workflow(ref_info, 'test_sample', fastq_input, {}, 1)
        self.assertGreater(wrapper.output.bam_file.size, 0)
        self.assertGreater(wrapper.output.vcf_unfiltered.size, 0)

    def test_variant_filtering(self) -> None:
        """
        Tests the variant filtering workflow.
        :return: None
        """
        wrapper = VariantFilteringWrapper(self.running_dir)
        wrapper.run_workflow(TestWorkflowsVariants.input_vcf, TestWorkflowsVariants.input_bam, {}, 1)
        self.assertGreater(Path(wrapper.output.vcf_filtered.path).stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
