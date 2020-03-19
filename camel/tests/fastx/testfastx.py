from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.fastx.fastqqualityfilter import FastqQualityFilter
from camel.app.tools.fastx.fastqqualitytrimmer import FastqQualityTrimmer


class TestFastx(CamelTestSuite):
    """
    Tests the fastx tools.
    """

    # Input files
    input_fastq = ToolIOFile(str(CamelTestSuite.get_test_file_dir('fastx') / 'ERR2019997-iontorrent-ds.fastq'))

    def test_quality_filter(self) -> None:
        """
        Tests the fastx quality filter.
        :return: None
        """
        q_filter = FastqQualityFilter(TestFastx.camel)
        q_filter.update_parameters(min_quality=30)
        q_filter.add_input_files({'FASTQ': [TestFastx.input_fastq]})
        q_filter.run(self.running_dir)
        self.assertGreater(TestFastx.input_fastq.size, q_filter.tool_outputs['FASTQ'][0].size)
        self.assertGreater(q_filter.informs['input_reads'], q_filter.informs['output_reads'])

    def test_quality_trimmer(self) -> None:
        """
        Tests the fastx quality trimmer.
        :return: None
        """
        q_trimmer = FastqQualityTrimmer(TestFastx.camel)
        q_trimmer.add_input_files({'FASTQ': [TestFastx.input_fastq]})
        q_trimmer.run(self.running_dir)
        self.assertGreater(TestFastx.input_fastq.size, q_trimmer.tool_outputs['FASTQ'][0].size)
        self.assertGreater(q_trimmer.informs['input_reads'], q_trimmer.informs['output_reads'])
