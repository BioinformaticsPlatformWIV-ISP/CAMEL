import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.spotyping.spotyping import SpoTyping


class TestSpoTyping(CamelTestSuite):
    """
    Tests the SpoTyping tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir()
    input_fastq_pe = [
        test_file_dir / 'pipelines' / 'Myco-DRR041783-ds_1.fastq.gz',
        test_file_dir / 'pipelines' / 'Myco-DRR041783-ds_2.fastq.gz'
    ]

    def test_spotyping(self) -> None:
        """
        Tests spotyping.
        :return: None
        """
        spotyping = SpoTyping(self.camel)
        spotyping.add_input_files({
            'FASTQ': [ToolIOFile(file_) for file_ in self.input_fastq_pe],
            'INFORMS_spoligo_param': [ToolIOValue({
                'min_strict': 3,
                'min_relaxed': 3,
                'downsample_factor': 'NA'
            })]
        })
        spotyping.run(self.running_dir)
        self.verify_output_files(spotyping, 'LOG')


if __name__ == '__main__':
    unittest.main()
