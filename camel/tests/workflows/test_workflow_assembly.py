from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.components.workflows.assemblywrapper import AssemblyWrapper
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
        test_file_dir / 'assembly' / 'reads_illumina_1.fastq',
        test_file_dir / 'assembly' / 'reads_illumina_2.fastq']
    fastq_se = test_file_dir / 'assembly' / 'reads_iontorrent.fastq'

    def test_assembly_illumina(self) -> None:
        """
        Tests the assembly workflow.
        :return: None
        """
        wrapper = AssemblyWrapper(self.running_dir)
        fastq_input = FastqInput('illumina', pe=[ToolIOFile(x) for x in TestWorkflowAssembly.fastq_pe])
        wrapper.run('test_sample', fastq_input, kmers='25,33')
        self.assertGreater(wrapper.output.fasta_contigs.stat().st_size, 0)

    def test_assembly_illumina_stats(self) -> None:
        """
        Tests the assembly workflow with standard forward / reverse reads input and with stats determined.
        :return: None
        """
        wrapper = AssemblyWrapper(self.running_dir)
        fastq_input = FastqInput('illumina', pe=[ToolIOFile(x) for x in TestWorkflowAssembly.fastq_pe])
        wrapper.run('test_sample', fastq_input, kmers='25,33', calc_qc_stats=True)
        self.assertGreater(wrapper.output.fasta_contigs.stat().st_size, 0)
        self.assertIsNotNone(wrapper.output.qc_stats)
        self.assertIn('depth', wrapper.output.qc_stats)
        self.assertIn('mapping', wrapper.output.qc_stats)

    def test_assembly_iontorrent(self) -> None:
        """
        Tests the assembly workflow.
        :return: None
        """
        wrapper = AssemblyWrapper(self.running_dir)
        fastq_input = FastqInput('iontorrent', se=[ToolIOFile(TestWorkflowAssembly.fastq_se)], is_pe=False)
        wrapper.run('test_sample', fastq_input, kmers='25,33', cov_cutoff=5, min_contig_length=500)
        self.assertGreater(wrapper.output.fasta_contigs.stat().st_size, 0)

    def test_assembly_iontorrent_stats(self) -> None:
        """
        Tests the assembly workflow with stats generation.
        :return: None
        """
        wrapper = AssemblyWrapper(self.running_dir)
        fastq_input = FastqInput('iontorrent', se=[ToolIOFile(TestWorkflowAssembly.fastq_se)], is_pe=False)
        wrapper.run('test_sample', fastq_input, kmers='25,33', cov_cutoff='auto', calc_qc_stats=True)
        self.assertGreater(wrapper.output.fasta_contigs.stat().st_size, 0)
        self.assertIsNotNone(wrapper.output.qc_stats)
        self.assertIn('depth', wrapper.output.qc_stats)
        self.assertIn('mapping', wrapper.output.qc_stats)
