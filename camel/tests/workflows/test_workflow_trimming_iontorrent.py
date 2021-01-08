from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.components.workflows.trimmingiontorrentwrapper import TrimmingIonTorrentWrapper


class TestWorkflowTrimmingIonTorrent(CamelTestSuite):
    """
    Tests the Illumina trimming workflow.
    """
    # Input files (Gene detection)
    test_file_dir = CamelTestSuite.get_test_file_dir('workflows')

    # Input files
    fastq_se = test_file_dir / 'trimming' / 'reads_iontorrent.fastq'
    fastq_se_gz = test_file_dir / 'trimming' / 'reads_iontorrent.fastq'

    def test_trimming_workflow_iontorrent(self) -> None:
        """
        Tests the read trimming workflow.
        :return: None
        """
        wrapper = TrimmingIonTorrentWrapper(self.running_dir)
        wrapper.run_workflow(TestWorkflowTrimmingIonTorrent.fastq_se)
        self.assertGreater(wrapper.output.trimmed_reads[0].size, 0)

    def test_trimming_workflow_iontorrent_gz(self) -> None:
        """
        Tests the read trimming workflow with gzipped input.
        :return: None
        """
        wrapper = TrimmingIonTorrentWrapper(self.running_dir)
        wrapper.run_workflow(TestWorkflowTrimmingIonTorrent.fastq_se_gz)
        self.assertGreater(wrapper.output.trimmed_reads[0].size, 0)
