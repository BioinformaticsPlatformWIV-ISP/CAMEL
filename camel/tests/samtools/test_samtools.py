import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.samtools.samtoolsdepth import SamtoolsDepth
from camel.app.tools.samtools.samtoolsdepthstatsanalyzer import SamtoolsDepthStatsAnalyzer
from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex
from camel.app.tools.samtools.samtoolsflagstat import SamtoolsFlagstat
from camel.app.tools.samtools.samtoolsindex import SamtoolsIndex
from camel.app.tools.samtools.samtoolsindexcram import SamtoolsIndexCram
from camel.app.tools.samtools.samtoolsmerge import SamtoolsMerge
from camel.app.tools.samtools.samtoolsmpileup import SamtoolsMPileup
from camel.app.tools.samtools.samtoolssort import SamtoolsSort
from camel.app.tools.samtools.samtoolsview import SamtoolsView


class TestSamtools(CamelTestSuite):
    """
    Tests the Samtools tool suite.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('samtools')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'reference.fasta')
    FILE_BAM = ToolIOFile(test_file_dir / 'aln.bam')
    FILE_BAM1 = ToolIOFile(test_file_dir / 'aln1.bam')
    FILE_BAM2 = ToolIOFile(test_file_dir / 'aln2.bam')
    FILE_CRAM = ToolIOFile(test_file_dir / 'aln.cram')
    FILE_TXTdepth = ToolIOFile(test_file_dir / 'aln_depth.txt')

    def test_samtools_depth(self) -> None:
        """
        Test SamtoolsDepth
        :return: None
        """
        samtools_depth = SamtoolsDepth(self.camel)
        samtools_depth.add_input_files({
            'BAM': [TestSamtools.FILE_BAM],
        })
        samtools_depth.run(self.running_dir)
        self.assertIn('TSV', samtools_depth.tool_outputs, "No TSV output generated.")
        output_file = samtools_depth.tool_outputs['TSV'][0].path
        self.assertTrue(output_file.exists(), "Output file does not exist.")
        self.assertGreater(output_file.stat().st_size, 0)
        self.assertIn('median_depth', samtools_depth.informs)

    def test_samtools_depthstatsanalyzer(self) -> None:
        """
        Tests SamtoolsDepthstatsanalyzer.
        :return: None
        """
        samtools_depthstatsanalyzer = SamtoolsDepthStatsAnalyzer(self.camel)

        samtools_depthstatsanalyzer.add_input_files({
            'TXT': [TestSamtools.FILE_TXTdepth],
            'FASTA_REF': [TestSamtools.FILE_FASTA_REF]
        })
        samtools_depthstatsanalyzer.run(self.running_dir)

        samtools_depthstatsanalyzer_expected_informs = ['median_coverage', 'coverage_iqr', 'coverage_cv', 'coverage_mad',
                                                        'coverage_std', 'base_coverage', 'segment_base_count', 'segment_gaps',
                                                        'segment_median_coverage', 'segment_coverage_mad', 'segment_coverage_cv',
                                                        'segment_coverage_iqr', 'segment_coverage_std', 'segment_base_coverage']

        for expected_inform in samtools_depthstatsanalyzer_expected_informs:
            self.assertIn(expected_inform, samtools_depthstatsanalyzer.informs, f"{expected_inform} not found in informs.")

    def test_samtools_fastaindex(self) -> None:
        """
        Tests SamtoolsFastaIndex.
        :return: None
        """
        samtools_fastaindex = SamtoolsFastaIndex(self.camel)

        samtools_fastaindex.add_input_files({
            'FASTA': [TestSamtools.FILE_FASTA_REF]
        })
        samtools_fastaindex.run(self.running_dir)

        self.assertIn('FASTA', samtools_fastaindex.tool_outputs, "No FASTA output generated.")
        output_file = samtools_fastaindex.tool_outputs['FASTA'][0].path
        self.assertTrue(output_file.exists(), "Output file does not exist.")
        self.assertGreater(output_file.stat().st_size, 0)

    def test_samtools_flagstat(self) -> None:
        """
        Tests SamtoolsFlagstat.
        :return: None
        """
        samtools_flagstat = SamtoolsFlagstat(self.camel)

        samtools_flagstat.add_input_files({
            'BAM': [TestSamtools.FILE_BAM]
        })
        samtools_flagstat.run(self.running_dir)

        self.assertIn('TXT', samtools_flagstat.tool_outputs, "No TXT output generated.")
        output_file = samtools_flagstat.tool_outputs['TXT'][0].path
        self.assertTrue(output_file.exists(), "Output file does not exist.")
        self.assertGreater(output_file.stat().st_size, 0)

        samtools_flagstat_expected_informs = ['total', 'secondary', 'supplementary', 'duplicates', 'mapped', 'paired',
                                              'read1', 'read2', 'properly_paired', 'singletons']

        for expected_inform in samtools_flagstat_expected_informs:
            self.assertIn(expected_inform, samtools_flagstat.informs, f"{expected_inform} not found in informs.")

    def test_samtools_index(self) -> None:
        """
        Tests SamtoolsIndex.
        :return: None
        """
        samtools_index = SamtoolsIndex(self.camel)
        samtools_index.add_input_files({
            'BAM': [TestSamtools.FILE_BAM],
        })
        samtools_index.run(self.running_dir)
        self.assertIn('BAM', samtools_index.tool_outputs, "No BAM output generated.")
        output_file = samtools_index.tool_outputs['BAM'][0].path
        self.assertTrue(output_file.exists(), "Output file does not exist.")
        self.assertGreater(output_file.stat().st_size, 0)

    def test_samtools_index_cram(self) -> None:
        """
        Tests SamtoolsIndexCram.
        :return: None
        """
        samtools_index_cram = SamtoolsIndexCram(self.camel)
        samtools_index_cram.add_input_files({
            'CRAM': [TestSamtools.FILE_CRAM],
            'FASTA_REF': [TestSamtools.FILE_FASTA_REF]
        })
        samtools_index_cram.run(self.running_dir)
        self.assertIn('CRAI', samtools_index_cram.tool_outputs, "No CRAI output generated.")
        output_file = samtools_index_cram.tool_outputs['CRAI'][0].path
        self.assertTrue(output_file.exists(), "Output file does not exist.")
        self.assertGreater(output_file.stat().st_size, 0)

    def test_samtools_merge(self) -> None:
        """
        Tests SamtoolsMerge.
        :return: None
        """
        samtools_merge = SamtoolsMerge(self.camel)
        samtools_merge.add_input_files({
            'BAM': [TestSamtools.FILE_BAM1, TestSamtools.FILE_BAM2],
        })
        samtools_merge.run(self.running_dir)
        self.assertIn('BAM', samtools_merge.tool_outputs, "No BAM output generated.")
        output_file = samtools_merge.tool_outputs['BAM'][0].path
        self.assertTrue(output_file.exists(), "Output file does not exist.")
        self.assertGreater(output_file.stat().st_size, 0)

    def test_samtools_mpileup(self) -> None:
        """
        Tests SamtoolsMPileup.
        :return: None
        """
        samtools_mpileup = SamtoolsMPileup(self.camel)
        samtools_mpileup.add_input_files({
            'BAM': [TestSamtools.FILE_BAM],
            'FASTA': [TestSamtools.FILE_FASTA_REF]
        })
        samtools_mpileup.run(self.running_dir)
        self.assertIn('PILEUP', samtools_mpileup.tool_outputs, "No PILEUP output generated.")
        output_file = samtools_mpileup.tool_outputs['PILEUP'][0].path
        self.assertTrue(output_file.exists(), "Output file does not exist.")
        self.assertGreater(output_file.stat().st_size, 0)

    def test_samtools_sort(self) -> None:
        """
        Tests SamtoolsSort.
        :return: None
        """
        samtools_sort = SamtoolsSort(self.camel)
        samtools_sort.add_input_files({
            'BAM': [TestSamtools.FILE_BAM],
        })
        samtools_sort.run(self.running_dir)
        self.assertIn('BAM', samtools_sort.tool_outputs, "No BAM output generated")
        output_file = samtools_sort.tool_outputs['BAM'][0].path
        self.assertTrue(output_file.exists(), "Output file does not exist.")
        self.assertGreater(output_file.stat().st_size, 0)

    def test_samtools_view(self) -> None:
        """
        Tests SamtoolsView.
        :return: None
        """
        samtools_view = SamtoolsView(self.camel)
        samtools_view.add_input_files({
            'BAM': [TestSamtools.FILE_BAM],
        })
        samtools_view.update_parameters(output_format = "SAM")
        samtools_view.run(self.running_dir)
        self.assertIn('SAM', samtools_view.tool_outputs, "No SAM output generated.")
        output_file = samtools_view.tool_outputs['SAM'][0].path
        self.assertTrue(output_file.exists(), "Output file does not exist.")
        self.assertGreater(output_file.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
