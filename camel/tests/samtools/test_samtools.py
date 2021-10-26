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
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'Homo_sapiens_assembly38_chr22.fasta')
    FILE_BAM = ToolIOFile(test_file_dir / 'NA12877_chr22_sub.bam')
    FILE_BAM2 = ToolIOFile(test_file_dir / 'NA12877_chr22_sub2.bam')
    FILE_CRAM = ToolIOFile(test_file_dir / 'NA12877_chr22_sub.cram')
    FILE_TXTdepth = ToolIOFile(test_file_dir / 'NA12877_chr22_sub_depth.txt')

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
        self.assertIn('TSV', samtools_depth.tool_outputs, "No TSV output generated")
        output_file = samtools_depth.tool_outputs['TSV'][0].path
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)
        self.assertIn('median_depth', samtools_depth.informs)

    def test_samtools_depthstatsanalyzer(self) -> None:
        """
        Test SamtoolsDepthstatsanalyzer
        :return: None
        """
        samtools_depthstatsanalyzer = SamtoolsDepthStatsAnalyzer(self.camel)

        samtools_depthstatsanalyzer.add_input_files({
            'TXT': [TestSamtools.FILE_TXTdepth],
            'FASTA_REF': [TestSamtools.FILE_FASTA_REF]
        })
        samtools_depthstatsanalyzer.run(self.running_dir)

        self.assertIn('median_coverage', samtools_depthstatsanalyzer.informs)
        self.assertIn('coverage_iqr', samtools_depthstatsanalyzer.informs)
        self.assertIn('coverage_cv', samtools_depthstatsanalyzer.informs)
        self.assertIn('coverage_mad', samtools_depthstatsanalyzer.informs)
        self.assertIn('coverage_std', samtools_depthstatsanalyzer.informs)
        self.assertIn('base_coverage', samtools_depthstatsanalyzer.informs)

        self.assertIn('segment_base_count', samtools_depthstatsanalyzer.informs)
        self.assertIn('segment_gaps', samtools_depthstatsanalyzer.informs)
        self.assertIn('segment_median_coverage', samtools_depthstatsanalyzer.informs)
        self.assertIn('segment_coverage_mad', samtools_depthstatsanalyzer.informs)
        self.assertIn('segment_coverage_cv', samtools_depthstatsanalyzer.informs)
        self.assertIn('segment_coverage_iqr', samtools_depthstatsanalyzer.informs)
        self.assertIn('segment_coverage_std', samtools_depthstatsanalyzer.informs)
        self.assertIn('segment_base_coverage', samtools_depthstatsanalyzer.informs)

    def test_samtools_fastaindex(self) -> None:
        """
        Test SamtoolsFastaIndex
        :return: None
        """
        samtools_fastaindex = SamtoolsFastaIndex(self.camel)

        samtools_fastaindex.add_input_files({
            'FASTA': [TestSamtools.FILE_FASTA_REF]
        })
        samtools_fastaindex.run(self.running_dir)

        self.assertIn('FASTA', samtools_fastaindex.tool_outputs, "No FASTA output generated")
        output_file = samtools_fastaindex.tool_outputs['FASTA'][0].path
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_samtools_flagstat(self) -> None:
        """
        Test SamtoolsFlagstat
        :return: None
        """
        samtools_flagstat = SamtoolsFlagstat(self.camel)

        samtools_flagstat.add_input_files({
            'BAM': [TestSamtools.FILE_BAM]
        })
        samtools_flagstat.run(self.running_dir)

        self.assertIn('TXT', samtools_flagstat.tool_outputs, "No TXT output generated")
        output_file = samtools_flagstat.tool_outputs['TXT'][0].path
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

        self.assertIn('total', samtools_flagstat.informs)
        self.assertIn('secondary', samtools_flagstat.informs)
        self.assertIn('supplementary', samtools_flagstat.informs)
        self.assertIn('duplicates', samtools_flagstat.informs)
        self.assertIn('mapped', samtools_flagstat.informs)
        self.assertIn('paired', samtools_flagstat.informs)
        self.assertIn('read1', samtools_flagstat.informs)
        self.assertIn('read2', samtools_flagstat.informs)
        self.assertIn('properly_paired', samtools_flagstat.informs)
        self.assertIn('singletons', samtools_flagstat.informs)

    def test_samtools_index(self) -> None:
        """
        Test SamtoolsIndex
        :return: None
        """
        samtools_index = SamtoolsIndex(self.camel)
        samtools_index.add_input_files({
            'BAM': [TestSamtools.FILE_BAM],
        })
        samtools_index.run(self.running_dir)
        self.assertIn('BAM', samtools_index.tool_outputs, "No BAM output generated")
        output_file = samtools_index.tool_outputs['BAM'][0].path
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_samtools_index_cram(self) -> None:
        """
        Test SamtoolsIndexCram
        :return: None
        """
        samtools_index_cram = SamtoolsIndexCram(self.camel)
        samtools_index_cram.add_input_files({
            'CRAM': [TestSamtools.FILE_CRAM],
            'FASTA_REF': [TestSamtools.FILE_FASTA_REF]
        })
        samtools_index_cram.run(self.running_dir)
        self.assertIn('CRAI', samtools_index_cram.tool_outputs, "No CRAI output generated")
        output_file = samtools_index_cram.tool_outputs['CRAI'][0].path
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_samtools_merge(self) -> None:
        """
        Test SamtoolsMerge
        :return: None
        """
        samtools_merge = SamtoolsMerge(self.camel)
        samtools_merge.add_input_files({
            'BAM': [TestSamtools.FILE_BAM, TestSamtools.FILE_BAM2],
        })
        samtools_merge.run(self.running_dir)
        self.assertIn('BAM', samtools_merge.tool_outputs, "No BAM output generated")
        output_file = samtools_merge.tool_outputs['BAM'][0].path
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_samtools_mpileup(self) -> None:
        """
        Test SamtoolsMPileup
        :return: None
        """
        samtools_mpileup = SamtoolsMPileup(self.camel)
        samtools_mpileup.add_input_files({
            'BAM': [TestSamtools.FILE_BAM],
            'FASTA': [TestSamtools.FILE_FASTA_REF]
        })
        samtools_mpileup.run(self.running_dir)
        self.assertIn('PILEUP', samtools_mpileup.tool_outputs, "No PILEUP output generated")
        output_file = samtools_mpileup.tool_outputs['PILEUP'][0].path
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_samtools_sort(self) -> None:
        """
        Test SamtoolsSort
        :return: None
        """
        samtools_sort = SamtoolsSort(self.camel)
        samtools_sort.add_input_files({
            'BAM': [TestSamtools.FILE_BAM],
        })
        samtools_sort.run(self.running_dir)
        self.assertIn('BAM', samtools_sort.tool_outputs, "No BAM output generated")
        output_file = samtools_sort.tool_outputs['BAM'][0].path
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_samtools_view(self) -> None:
        """
        Test SamtoolsView
        :return: None
        """
        samtools_view = SamtoolsView(self.camel)
        samtools_view.add_input_files({
            'BAM': [TestSamtools.FILE_BAM],
        })
        samtools_view.update_parameters(output_format = "SAM")
        samtools_view.run(self.running_dir)
        self.assertIn('SAM', samtools_view.tool_outputs, "No SAM output generated")
        output_file = samtools_view.tool_outputs['SAM'][0].path
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
