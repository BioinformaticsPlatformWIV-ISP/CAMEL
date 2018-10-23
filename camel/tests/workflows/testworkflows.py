import os
import unittest

import tempfile

from camel.app.camel import Camel
from camel.app.components.workflows.assemblywrapper import AssemblyWrapper
from camel.app.components.workflows.genedetectionwrapper import GeneDetectionWrapper
from camel.app.components.workflows.readtrimmingwrapper import ReadTrimmingWrapper
from camel.app.io.tooliofile import ToolIOFile


class TestWorkflows(unittest.TestCase):
    """
    Tests the Snakemake workflows.
    """
    camel = Camel()
    running_dir = None

    # Input files
    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'])
    input_fasta = ToolIOFile(os.path.join(test_file_dir, 'workflows', 'NC_002695.1.fasta'))
    input_reads_raw = [ToolIOFile(os.path.join(test_file_dir, 'workflows', 'ecoli_1.fastq')),
                       ToolIOFile(os.path.join(test_file_dir, 'workflows', 'ecoli_2.fastq'))]
    input_reads_trim = {
        'PE': [ToolIOFile(os.path.join(test_file_dir, 'workflows', 'ecoli_trim_1.fastq')),
               ToolIOFile(os.path.join(test_file_dir, 'workflows', 'ecoli_trim_2.fastq'))],
        'FWD': [ToolIOFile(os.path.join(test_file_dir, 'workflows', 'ecoli_trim_fwd_only.fastq'))],
        'REV': [ToolIOFile(os.path.join(test_file_dir, 'workflows', 'ecoli_trim_rev_only.fastq'))]
    }
    input_gene_detection_db = os.path.join(test_file_dir, 'gene_detection', 'db')

    def setUp(self) -> None:
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(prefix='camel_', dir=TestWorkflows.camel.config['temp_dir'])

    def test_trimming_workflow(self) -> None:
        """
        Tests the read trimming workflow.
        :return: None
        """
        wrapper = ReadTrimmingWrapper(self.running_dir)
        wrapper.run_workflow([x.path for x in TestWorkflows.input_reads_raw])
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
        wrapper.run_workflow('test_sample', TestWorkflows.input_reads_raw, [], [], kmers='25,29,33')
        self.assertGreater(wrapper.output.fasta_contigs.size, 0)

    def test_assembly_workflow_trimmed(self) -> None:
        """
        Tests the assembly workflow on trimmed reads.
        :return: None
        """
        wrapper = AssemblyWrapper(self.running_dir)
        wrapper.run_workflow('test_sample', TestWorkflows.input_reads_trim['PE'],
                             TestWorkflows.input_reads_trim['FWD'],
                             TestWorkflows.input_reads_trim['REV'],
                             kmers='55')
        self.assertGreater(wrapper.output.fasta_contigs.size, 0)

    def test_gene_detection_workflow_blast(self) -> None:
        """
        Tests the gene detection workflow using BLAST.
        :return: None
        """
        wrapper = GeneDetectionWrapper(self.running_dir)
        db_data = {
            'path': TestWorkflows.input_gene_detection_db,
            'min_percent_identity': 80,
            'min_coverage': 80
        }
        wrapper.run_workflow_blast(TestWorkflows.input_fasta.path, 'test_sample', db_data)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)

    def test_gene_detection_workflow_srst2(self) -> None:
        """
        Tests the gene detection workflow using SRST2.
        :return: None
        """
        wrapper = GeneDetectionWrapper(self.running_dir)
        db_data = {
            'path': TestWorkflows.input_gene_detection_db,
            'min_percent_identity': 80,
            'min_coverage': 80
        }
        wrapper.run_workflow_srst2([x.path for x in TestWorkflows.input_reads_raw], 'test_sample', db_data)
        self.assertGreater(len(wrapper.output.report_section.to_html()), 0)
