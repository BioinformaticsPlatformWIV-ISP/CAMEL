import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.ncbihumanreadscrubber.ncbihumanreadscrubber import NcbiHumanReadScrubber
from camel.scripts.ncbihumanreadscrubber.mainncbihumanreadscrubber import MainNcbiHumanReadScrubber

class TestNcbiHumanReadScrubber(CamelTestSuite):
    """
    Tests the HRRT tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('ncbihumanreadscrubber')

    # def test_scrubber(self) -> None:
    #     """
    #     Tests the scrubber on a not gzipped single end fastq file.
    #     :return: None
    #     """
    #     scrubber = NcbiHumanReadScrubber(self.camel)
    #     scrubber.add_input_files({'FASTQ_SINGLE_GUNZIP':
    #                               [ToolIOFile(TestNcbiHumanReadScrubber.test_file_dir / 'reads_illumina_1.fastq')]})
    #     scrubber.update_parameters(interleaved='false', outputfile=self.running_dir / 'test_scrubber_output.fastq')
    #     scrubber.run(self.running_dir)
    #     self.verify_output_files(scrubber, 'FASTQ_SCRUBBED', 1)

    def test_scrubbing_pipeline(self) -> None:
        """
        Tests the NCBI human read scrubbing standalone pipeline.
        :return: None
        """
        path_report_out = self.running_dir / 'out' / 'report.html'
        args = [
            '--fastq-pe', str(TestNcbiHumanReadScrubber.test_file_dir / 'reads_illumina_1.fastq'), str(TestNcbiHumanReadScrubber.test_file_dir / 'reads_illumina_2.fastq'),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent),
            '--working-dir', str(self.running_dir),
            '--output-tsv', "None",
            '--read-type', 'illumina'
        ]
        main = MainNcbiHumanReadScrubber(args)
        main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
