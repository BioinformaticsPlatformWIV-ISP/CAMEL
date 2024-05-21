import unittest

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
    input_fasta = test_file_dir / 'pipelines' / 'Myco-DRR041783-ds.fasta'

    def test_spotyping(self) -> None:
        """
        Tests the SpoTyping tool.
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
        self.assertIn('VAL_type_binary', spotyping.tool_outputs)
        self.verify_output_files(spotyping, 'LOG')

    def test_spotyping_swift_off(self) -> None:
        """
        Tests the SpoTyping tool with the 'swift' parameter disabled.
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
        spotyping.update_parameters(swift='off')
        spotyping.run(self.running_dir)
        self.assertIn('VAL_type_binary', spotyping.tool_outputs)
        self.verify_output_files(spotyping, 'LOG')

    def test_spotyping_fasta(self) -> None:
        """
        Tests the SpoTyping tool with FASTA input.
        :return: None
        """
        spotyping = SpoTyping(self.camel)
        spotyping.add_input_files({
            'FASTA': [ToolIOFile(self.input_fasta)]
        })
        spotyping.update_parameters(fasta=None)
        spotyping.run(self.running_dir)
        self.assertIn('VAL_type_binary', spotyping.tool_outputs)
        self.verify_output_files(spotyping, 'LOG')


if __name__ == '__main__':
    unittest.main()
