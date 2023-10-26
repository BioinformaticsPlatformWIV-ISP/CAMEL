import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.ncbihumanreadscrubber.ncbihumanreadscrubber import NcbiHumanReadScrubber


class TestHRRT(CamelTestSuite):
    """
    Tests the HRRT tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('ncbihumanreadscrubber')

    def test_scrubber(self) -> None:
        """
        Tests the scrubber on a not gzipped single end fastq file.
        :return: None
        """
        scrubber = NcbiHumanReadScrubber(self.camel)
        scrubber.add_input_files({'FASTQ_INTERLEAVED_GUNZIPPED':
                                  [ToolIOFile(TestHRRT.test_file_dir / 'reads_illumina_1.fastq')]})
        scrubber.update_parameters(interleaved='false', outputfile=self.running_dir / 'test_scrubber_output.fastq')
        scrubber.run(self.running_dir)
        self.verify_output_files(scrubber, 'FASTQ_SCRUBBED', 1)


if __name__ == '__main__':
    unittest.main()
