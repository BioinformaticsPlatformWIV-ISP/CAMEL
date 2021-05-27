from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.components.workflows.trimmingilluminawrapper import TrimmingIlluminaWrapper


class TestWorkflowTrimmingIllumina(CamelTestSuite):
    """
    Tests the Illumina trimming workflow.
    """
    # Input files (Gene detection)
    test_file_dir = CamelTestSuite.get_test_file_dir('workflows')

    # Input files
    fastq_pe = [test_file_dir / 'trimming' / 'reads_1.fastq', test_file_dir / 'trimming' / 'reads_2.fastq']
    fastq_pe_gz = [test_file_dir / 'trimming' / 'reads_1.fastq.gz', test_file_dir / 'trimming' / 'reads_2.fastq.gz']

    def test_trimming_workflow_illumina(self) -> None:
        """
        Tests the read trimming workflow.
        :return: None
        """
        wrapper = TrimmingIlluminaWrapper(self.running_dir)
        wrapper.run_workflow(TestWorkflowTrimmingIllumina.fastq_pe)
        self.assertGreater(wrapper.output.trimmed_reads_pe[0].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_pe[1].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_se_fwd[0].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_se_rev[0].size, 0)

    def test_trimming_workflow_illumina_truseq(self) -> None:
        """
        Tests the read trimming workflow with TruSeq adapters.
        :return: None
        """
        wrapper = TrimmingIlluminaWrapper(self.running_dir)
        wrapper.run_workflow(TestWorkflowTrimmingIllumina.fastq_pe, 'TruSeq3')
        self.assertGreater(wrapper.output.trimmed_reads_pe[0].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_pe[1].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_se_fwd[0].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_se_rev[0].size, 0)
        self.assertIn('TruSeq3', wrapper.output.informs_trimmomatic['_command'])

    def test_trimming_workflow_illumina_gz(self) -> None:
        """
        Tests the read trimming workflow with gzipped input.
        :return: None
        """
        wrapper = TrimmingIlluminaWrapper(self.running_dir)
        wrapper.run_workflow(TestWorkflowTrimmingIllumina.fastq_pe_gz)
        self.assertGreater(wrapper.output.trimmed_reads_pe[0].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_pe[1].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_se_fwd[0].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_se_rev[0].size, 0)
