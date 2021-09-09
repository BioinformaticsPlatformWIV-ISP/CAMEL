import unittest
from pathlib import Path

from camel.app.camel import Camel
from camel.app.components.pipelines import pipeutils
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.bowtie2.bowtie2map import Bowtie2Map
from camel.app.tools.samtools.samtoolsflagstat import SamtoolsFlagstat
from camel.app.tools.samtools.samtoolssort import SamtoolsSort
from camel.app.tools.samtools.samtoolsview import SamtoolsView


class TestPipedTools(CamelTestSuite):
    """
    Tests the piping of tools.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('read_mapping')
    input_reference = test_file_dir / 'reference.fasta'
    input_fq_pe = [
        test_file_dir / 'reads_1.fastq',
        test_file_dir / 'reads_2.fastq'
    ]

    def test_pipeline_bt2_sam_to_bam(self) -> None:
        """
        Tests the tree construction.
        :return: None
        """
        # Initialize tools
        bowtie2 = Bowtie2Map(Camel.get_instance())
        bowtie2.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in TestPipedTools.input_fq_pe],
            'INDEX_GENOME_PREFIX': [ToolIOValue(TestPipedTools.input_reference)]
        })
        samtools_view = SamtoolsView(Camel.get_instance())

        # Run as pipe
        pipeutils.run_as_pipe([bowtie2, samtools_view], self.running_dir)

        # Assert output file is generated
        self.assertIn('BAM', samtools_view.tool_outputs)
        self.assertTrue(Path(samtools_view.tool_outputs['BAM'][0].path))
        self.assertGreater(Path(samtools_view.tool_outputs['BAM'][0].path).stat().st_size, 0)

        # Assert Bowtie2 informs are parsed
        self.assertIn('stats_map_rate', bowtie2.informs)

    def test_pipeline_bt2_sam_to_bam_invalid_input(self) -> None:
        """
        Tests if the pipeline raises an error if the input of t he workflow is invalid.
        :return: None
        """
        # Initialize tools
        bowtie2 = Bowtie2Map(Camel.get_instance())
        bowtie2.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in ('non_existing_file')],
            'INDEX_GENOME_PREFIX': [ToolIOValue(TestPipedTools.input_reference)]
        })
        samtools_view = SamtoolsView(Camel.get_instance())

        # Run as pipe
        with self.assertRaises(BaseException):
            pipeutils.run_as_pipe([bowtie2, samtools_view], self.running_dir)

    def test_pipeline_bt2_sam_to_sorted_bam(self) -> None:
        """
        Tests the tree construction.

        :return: None
        """
        # Initialize tools
        bowtie2 = Bowtie2Map(Camel.get_instance())
        bowtie2.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in TestPipedTools.input_fq_pe],
            'INDEX_GENOME_PREFIX': [ToolIOValue(TestPipedTools.input_reference)]
        })
        samtools_view = SamtoolsView(Camel.get_instance())
        samtools_sort = SamtoolsSort(Camel.get_instance())

        # Run as pipe
        pipeutils.run_as_pipe([bowtie2, samtools_view, samtools_sort], self.running_dir)

        # Assert output file is generated
        self.assertIn('BAM', samtools_sort.tool_outputs)
        self.assertTrue(Path(samtools_sort.tool_outputs['BAM'][0].path))
        self.assertGreater(Path(samtools_sort.tool_outputs['BAM'][0].path).stat().st_size, 0)

        # Assert that intermediate step do not generate output
        self.assertNotIn('BAM', samtools_view.tool_outputs)

        # Assert Bowtie2 informs are parsed
        self.assertIn('stats_map_rate', bowtie2.informs)

    def test_pipeline_bt2_sam_to_bam_to_flagstat(self) -> None:
        """
        Tests the tree construction.
        :return: None
        """
        # Initialize tools
        bowtie2 = Bowtie2Map(Camel.get_instance())
        bowtie2.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in TestPipedTools.input_fq_pe],
            'INDEX_GENOME_PREFIX': [ToolIOValue(TestPipedTools.input_reference)]
        })
        samtools_view = SamtoolsView(Camel.get_instance())
        samtools_sort = SamtoolsSort(Camel.get_instance())
        samtools_flagstat = SamtoolsFlagstat(Camel.get_instance())

        # Run as pipe
        pipeutils.run_as_pipe([bowtie2, samtools_view, samtools_sort, samtools_flagstat], self.running_dir)

        # Assert samtools sort does not generate output
        self.assertIn('TXT', samtools_flagstat.tool_outputs)
        self.assertGreater(Path(samtools_flagstat.tool_outputs['TXT'][0].path).stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
