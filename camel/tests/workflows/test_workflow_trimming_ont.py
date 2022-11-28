import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.components.workflows.trimmingontwrapper import TrimmingONTWrapper


class TestWorkflowTrimmingONT(CamelTestSuite):
    """
    Tests the ONT trimming workflow.
    """
    # Retrieving the filetest directory
    test_file_dir = CamelTestSuite.get_test_file_dir('workflows')

    # Input files
    fastq_se = test_file_dir / 'trimming' / 'subsampled_ont_13_2500.fastq'
    fastq_se_gz = test_file_dir / 'trimming' / 'subsampled_ont_13_2500.fastq.gz'

    def test_trimming_workflow_ont(self) -> None:
        """
        Tests the read trimming workflow.
        :return: None
        """
        wrapper = TrimmingONTWrapper(self.running_dir)
        wrapper.run_workflow(TestWorkflowTrimmingONT.fastq_se)
        self.assertGreater(wrapper.output.trimmed_reads[0].size, 0)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)

    def test_trimming_workflow_ont_gz(self) -> None:
        """
        Tests the read trimming workflow with gzipped input.
        :return: None
        """
        wrapper = TrimmingONTWrapper(self.running_dir)
        wrapper.run_workflow(TestWorkflowTrimmingONT.fastq_se_gz)
        self.assertGreater(wrapper.output.trimmed_reads[0].size, 0)


if __name__ == '__main__':
    unittest.main()

