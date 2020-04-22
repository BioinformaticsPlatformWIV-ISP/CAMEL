import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.components.workflows.assemblywrapper import AssemblyWrapper
from camel.app.components.workflows.genedetectionwrapper import GeneDetectionWrapper
from camel.app.components.workflows.readtrimmingwrapper import ReadTrimmingWrapper
from camel.app.components.workflows.sequencetypingwrapper import SequenceTypingWrapper, SequenceTypingInput
from camel.app.io.tooliofile import ToolIOFile
from camel.tests import longRunningTest


class TestWorkflows(CamelTestSuite):
    """
    Tests the Snakemake workflows.
    """
    # Input files (Gene detection)
    test_file_dir = CamelTestSuite.get_test_file_dir('workflows')
    input_gene_fasta = test_file_dir / 'NC_002695.1.fasta'
    input_gene_reads_raw = [
        test_file_dir / 'ecoli_1.fastq',
        test_file_dir / 'ecoli_2.fastq']
    input_gene_reads_trim = {
        'PE': [test_file_dir / 'ecoli_trim_1.fastq', test_file_dir / 'ecoli_trim_2.fastq'],
        'FWD': [test_file_dir / 'ecoli_trim_fwd_only.fastq'],
        'REV': [test_file_dir / 'ecoli_trim_rev_only.fastq']
    }
    input_gene_db = CamelTestSuite.get_test_file_dir('gene_detection') / 'db'

    # Input files (typing)
    test_file_dir_typing = CamelTestSuite.get_test_file_dir('typing')
    input_typing_fasta = test_file_dir_typing / 'neisseria_mc58.fasta'
    input_typing_db = test_file_dir_typing / 'scheme_mlst_neisseria'
    input_typing_db_protein = test_file_dir_typing / 'scheme_pora_neisseria'
    input_typing_db_mixed = test_file_dir_typing / 'scheme_fhbp_neisseria'
    input_typing_reads = [
        test_file_dir_typing / 'S15BD05018_S58_L001_1.fastq',
        test_file_dir_typing / 'S15BD05018_S58_L001_1.fastq'
    ]

    def test_trimming_workflow(self) -> None:
        """
        Tests the read trimming workflow.
        :return: None
        """
        wrapper = ReadTrimmingWrapper(self.running_dir)
        wrapper.run_workflow([str(f) for f in TestWorkflows.input_gene_reads_raw])
        self.assertGreater(wrapper.output.trimmed_reads_pe[0].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_pe[1].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_se_fwd[0].size, 0)
        self.assertGreater(wrapper.output.trimmed_reads_se_rev[0].size, 0)

    def test_assembly_workflow_pe(self) -> None:
        """
        Tests the assembly workflow with standard forward / reverse reads input.
        :return: None
        """
        wrapper = AssemblyWrapper(self.running_dir)
        reads = [ToolIOFile(str(x)) for x in TestWorkflows.input_gene_reads_raw]
        wrapper.run_workflow('test_sample', reads, [], [], kmers='25,29,33', cov_cutoff=11)
        self.assertGreater(wrapper.output.fasta_contigs.size, 0)

    def test_assembly_workflow_trimmed(self) -> None:
        """
        Tests the assembly workflow on trimmed reads.
        :return: None
        """
        wrapper = AssemblyWrapper(self.running_dir)
        wrapper.run_workflow(
            'test_sample',
            [ToolIOFile(str(x)) for x in TestWorkflows.input_gene_reads_trim['PE']],
            [ToolIOFile(str(TestWorkflows.input_gene_reads_trim['FWD'][0]))],
            [ToolIOFile(str(TestWorkflows.input_gene_reads_trim['REV'][0]))],
            kmers='55')
        self.assertGreater(wrapper.output.fasta_contigs.size, 0)

    def test_assembly_workflow_stats(self) -> None:
        """
        Tests the assembly workflow with standard forward / reverse reads input and with stats determined.
        :return: None
        """
        wrapper = AssemblyWrapper(self.running_dir)
        reads = [ToolIOFile(str(x)) for x in TestWorkflows.input_gene_reads_raw]
        wrapper.run_workflow('test_sample', reads, [], [], kmers='25,29,33', cov_cutoff=11, calculate_qc_stats=True)
        self.assertGreater(wrapper.output.fasta_contigs.size, 0)
        self.assertIsNotNone(wrapper.output.qc_stats)
        self.assertIn('depth', wrapper.output.qc_stats)
        self.assertIn('mapping', wrapper.output.qc_stats)

    def test_gene_detection_workflow_blast(self) -> None:
        """
        Tests the gene detection workflow using BLAST.
        :return: None
        """
        wrapper = GeneDetectionWrapper(self.running_dir)
        db_data = {
            'path': str(TestWorkflows.input_gene_db),
            'min_percent_identity': 80,
            'min_coverage': 80
        }
        wrapper.run_workflow_blast(str(TestWorkflows.input_gene_fasta), 'test_sample', db_data)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)

    def test_gene_detection_workflow_srst2(self) -> None:
        """
        Tests the gene detection workflow using SRST2.
        :return: None
        """
        wrapper = GeneDetectionWrapper(self.running_dir)
        db_data = {
            'path': str(TestWorkflows.input_gene_db),
            'min_percent_identity': 80,
            'min_coverage': 80
        }
        wrapper.run_workflow_srst2([str(x) for x in TestWorkflows.input_gene_reads_raw], 'test_sample', db_data)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)

    def test_gene_detection_workflow_kma(self) -> None:
        """
        Tests the gene detection workflow using SRST2.
        :return: None
        """
        wrapper = GeneDetectionWrapper(self.running_dir)
        db_data = {'path': str(TestWorkflows.input_gene_db)}
        wrapper.run_workflow_kma([str(x) for x in TestWorkflows.input_gene_reads_raw], 'test_sample', db_data)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)

    def test_typing_workflow_blast(self) -> None:
        """
        Tests the sequence typing workflow using BLAST+.
        :return: None
        """
        wrapper = SequenceTypingWrapper(self.running_dir)
        workflow_input = SequenceTypingInput(
            sample_name='test_sample',
            db_path=str(TestWorkflows.input_typing_db),
            fasta=ToolIOFile(str(TestWorkflows.input_typing_fasta))
        )
        wrapper.run_workflow_blast(workflow_input, 8)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)

    def test_typing_workflow_blast_protein(self) -> None:
        """
        Tests the sequence typing workflow using BLAST+ with a protein scheme.
        :return: None
        """
        wrapper = SequenceTypingWrapper(self.running_dir)
        workflow_input = SequenceTypingInput(
            sample_name='test_sample',
            db_key='pora',
            db_path=str(TestWorkflows.input_typing_db_protein),
            fasta=ToolIOFile(str(TestWorkflows.input_typing_fasta))
        )
        wrapper.run_workflow_blast(workflow_input, 8)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)

    @longRunningTest()
    def test_typing_workflow_blast_mixed(self) -> None:
        """
        Tests the sequence typing workflow using BLAST+ with a mixed scheme.
        :return: None
        """
        wrapper = SequenceTypingWrapper(self.running_dir)
        workflow_input = SequenceTypingInput(
            sample_name='test_sample',
            db_key='fhbp',
            db_path=str(TestWorkflows.input_typing_db_mixed),
            fasta=ToolIOFile(str(TestWorkflows.input_typing_fasta))
        )
        wrapper.run_workflow_blast(workflow_input, 8)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)

    @longRunningTest()
    def test_typing_workflow_srst2(self) -> None:
        """
        Tests the sequence typing workflow using SRST2.
        :return: None
        """
        wrapper = SequenceTypingWrapper(self.running_dir)
        workflow_input = SequenceTypingInput(
            sample_name='test_sample',
            db_path=str(TestWorkflows.input_typing_db),
            fastq_pe=[ToolIOFile(str(x)) for x in TestWorkflows.input_typing_reads]
        )
        wrapper.run_workflow_srst2(workflow_input, {'max_unaligned_overlap': 100}, 8)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)

    @longRunningTest()
    def test_typing_workflow_srst2_mixed(self) -> None:
        """
        Tests the sequence typing workflow using SRST2 with a mixed scheme (protein and DNA loci).
        :return: None
        """
        wrapper = SequenceTypingWrapper(self.running_dir)
        workflow_input = SequenceTypingInput(
            sample_name='test_sample',
            db_path=str(TestWorkflows.input_typing_db),
            fastq_pe=[ToolIOFile(str(x)) for x in TestWorkflows.input_typing_reads],
            fasta=ToolIOFile(str(TestWorkflows.input_typing_fasta))
        )
        wrapper.run_workflow_srst2(workflow_input, threads=8)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
