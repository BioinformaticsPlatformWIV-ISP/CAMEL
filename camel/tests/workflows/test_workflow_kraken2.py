from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.components.workflows.assemblywrapper import AssemblyWrapper
from camel.app.components.workflows.kraken2wrapper import Kraken2Wrapper
from camel.app.components.workflows.utils.fastqinput import FastqInput
from camel.app.io.tooliofile import ToolIOFile


class TestWorkflowAssembly(CamelTestSuite):
    """
    Tests the Illumina trimming workflow.
    """
    # Input files (Gene detection)
    test_file_dir = CamelTestSuite.get_test_file_dir('workflows')

    # Input files
    fastq_pe = [
        test_file_dir / 'kraken2' / 'reads_illumina_1.fastq',
        test_file_dir / 'kraken2' / 'reads_illumina_2.fastq']
    fastq_se = test_file_dir / 'kraken2' / 'reads_iontorrent.fastq'

    def test_kraken2_illumina_paired_end(self) -> None:
        """
        Tests the KRAKEN2 workflow on Illumina PE data.
        :return: None
        """
        wrapper = Kraken2Wrapper(self.running_dir)
        fastq_input = FastqInput('illumina', pe=[ToolIOFile(x) for x in TestWorkflowAssembly.fastq_pe])
        expected_species = 'Escherichia coli'
        wrapper.run_workflow('test_sample', fastq_input, expected_species)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)
        self.assertGreater(wrapper.output.tsv_summary.stat().st_size, 0)

    def test_kraken2_iontorrent_single_end(self) -> None:
        """
        Tests the KRAKEN2 workflow on IonTorrent SE data.
        :return: None
        """
        wrapper = Kraken2Wrapper(self.running_dir)
        fastq_input = FastqInput('iontorrent', se=[ToolIOFile(TestWorkflowAssembly.fastq_se)], is_pe=False)
        expected_species = 'Escherichia coli'
        wrapper.run_workflow('test_sample', fastq_input, expected_species)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)
        self.assertGreater(wrapper.output.tsv_summary.stat().st_size, 0)

    def test_kraken2_iontorrent_single_end(self) -> None:
        """
        Tests the KRAKEN2 workflow on IonTorrent SE data.
        :return: None
        """
        wrapper = Kraken2Wrapper(self.running_dir)
        fastq_input = FastqInput('iontorrent', se=[ToolIOFile(TestWorkflowAssembly.fastq_se)], is_pe=False)
        expected_species = 'Escherichia coli'
        wrapper.run_workflow('test_sample', fastq_input, expected_species)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)
        self.assertGreater(wrapper.output.tsv_summary.stat().st_size, 0)

    def test_kraken2_illumina_paired_end_contaminated(self) -> None:
        """
        Tests the KRAKEN2 workflow on Illumina PE data, with reads other than the expected species.
        :return: None
        """
        wrapper = Kraken2Wrapper(self.running_dir)
        fastq_input = FastqInput('illumina', pe=[ToolIOFile(x) for x in TestWorkflowAssembly.fastq_pe])
        expected_species = 'Listeria monocytogenes'
        wrapper.run_workflow('test_sample', fastq_input, expected_species)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)
        self.assertGreater(wrapper.output.tsv_summary.stat().st_size, 0)
        self.assertGreater(len(wrapper.output.informs['contaminants_warn']), 0)

    def test_kraken2_illumina_paired_end_genus(self) -> None:
        """
        Tests the KRAKEN2 workflow on Illumina PE data.
        :return: None
        """
        wrapper = Kraken2Wrapper(self.running_dir)
        fastq_input = FastqInput('illumina', pe=[ToolIOFile(x) for x in TestWorkflowAssembly.fastq_pe])
        expected_species = 'Escherichia'
        wrapper.run_workflow('test_sample', fastq_input, expected_species, level_of_depth='G')
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)
        self.assertGreater(wrapper.output.tsv_summary.stat().st_size, 0)