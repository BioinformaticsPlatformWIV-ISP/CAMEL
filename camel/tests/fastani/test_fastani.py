import unittest
from pathlib import Path

from camel.app.camel import Camel
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.fastani.fastani import FastANI
from camel.app.tools.fastani.fastanireporter import FastANIReporter


class TestFastANI(CamelTestSuite):
    """
    Tests the FastANI tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('fastani')
    FILE_FASTA_R = ToolIOFile(test_file_dir / 'escherichia.fasta')
    FILE_FASTA_Q = ToolIOFile(test_file_dir / 'shigella.fasta')

    def test_fastani(self) -> None:
        """
        Tests the FastANI tool with default options.
        :return: None
        """
        fastani = FastANI(self.camel)
        fastani.add_input_files({'FASTA_Q': [TestFastANI.FILE_FASTA_Q], 'FASTA_R': [TestFastANI.FILE_FASTA_R]})
        fastani.run(self.running_dir)
        self.verify_output_files(fastani, 'TSV')

    def test_fastani_reporter(self) -> None: # Unfinished yet
        """
        Tests the FastANI reporter class.
        :return: None
        """
        # Run FastANI
        fastani = FastANI(self.camel)
        fastani.add_input_files({'FASTA': [TestFastANI.FILE_FASTA]})
        fastani.run(self.running_dir)
        self.verify_output_files(fastani, 'TSV')

        # Run the FastANI reporter
        reporter = FastANIReporter(self.camel)
        reporter.add_input_files({'TSV': fastani.tool_outputs['TSV']})
        reporter.add_input_informs({'fastani': fastani.informs})
        reporter.run(self.running_dir)
        output_section = reporter.tool_outputs['HTML'][0].value
        self.assertGreater(len(output_section.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
