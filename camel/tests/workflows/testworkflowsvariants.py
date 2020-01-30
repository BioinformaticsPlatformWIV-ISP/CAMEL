import unittest
from pathlib import Path

import tempfile

from camel.app.camel import Camel
from camel.app.components.workflows.variantcallingwrapper import VariantCallingWrapper, VariantCallingInput
from camel.app.components.workflows.variantfilteringwrapper import VariantFilteringWrapper
from camel.app.io.tooliofile import ToolIOFile


class TestWorkflowsVariants(unittest.TestCase):
    """
    Tests the Snakemake variant calling related workflows.
    """
    camel = Camel.get_instance()
    running_dir = None

    test_file_dir = Path(camel.config['testing']['testfiles_dir'])

    # Variant calling
    input_fastq_pe = [
        test_file_dir / 'variant_calling' / 'toy' / 'r_1.fq',
        test_file_dir / 'variant_calling' / 'toy' / 'r_2.fq']
    input_fasta_ref = test_file_dir / 'variant_calling' / 'toy' / 'toy.fasta'

    # Variant filtering
    input_vcf = test_file_dir / 'variant_calling' / 'unfiltered_variants-myco.vcf.gz'
    input_bam = test_file_dir / 'variant_calling' / 'alignment.bam'

    def setUp(self) -> None:
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = Path(tempfile.mkdtemp(prefix='camel_', dir=TestWorkflowsVariants.camel.config['temp_dir']))

    def test_variant_calling(self) -> None:
        """
        Tests the variant calling workflow.
        :return: None
        """
        wrapper = VariantCallingWrapper(self.running_dir)
        in_files = VariantCallingInput(pe_reads=[ToolIOFile(x) for x in TestWorkflowsVariants.input_fastq_pe])
        ref_info = {
            'name': TestWorkflowsVariants.input_fasta_ref.name,
            'path': str(TestWorkflowsVariants.input_fasta_ref)
        }
        wrapper.run_workflow(ref_info, 'test_sample', in_files, {}, 1)
        self.assertGreater(wrapper.output.bam_file.size, 0)
        self.assertGreater(wrapper.output.vcf_unfiltered.size, 0)

    def test_variant_filtering(self) -> None:
        """
        Tests the variant filtering workflow.
        :return: None
        """
        wrapper = VariantFilteringWrapper(self.running_dir)
        wrapper.run_workflow(str(TestWorkflowsVariants.input_vcf), TestWorkflowsVariants.input_bam, {}, 1)
        self.assertGreater(Path(wrapper.output.vcf_filtered.path).stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
