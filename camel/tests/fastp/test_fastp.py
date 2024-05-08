import unittest

from camel.app.camel import Camel
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.fastp.fastp import Fastp


class TestFastp(CamelTestSuite):
    """
    Tests for the fastp tool class.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('fastp')
    INPUT_FASTQ_PE = [
        ToolIOFile(test_file_dir / 'reads_illumina_1.fastq.gz'),
        ToolIOFile(test_file_dir / 'reads_illumina_2.fastq.gz')]

    def test_seqkit_fastp_pe(self) -> None:
        """
        Tests fastp with PE input.
        """
        fastp = Fastp(self.camel)
        fastp.add_input_files({'FASTQ': TestFastp.INPUT_FASTQ_PE})
        fastp.update_parameters(
            output_name='reads_illumina',
            # Adapter trimming
            detect_adapter_for_pe=True,
            # Leading
            cut_front=True,
            cut_front_window_size=1,
            cut_front_mean_quality=10,
            # Trailing
            cut_tail=True,
            cut_tail_window_size=1,
            cut_tail_mean_quality=10,
            # Sliding window
            cut_right=True,
            cut_right_window_size=4,
            cut_right_mean_quality=20,
            # Minimum length
            length_required=40
        )
        fastp.run(self.running_dir)
        self.verify_output_files(fastp, 'FASTQ_PE', nb_files=2)
        self.verify_output_files(fastp, 'JSON')
        self.verify_output_files(fastp, 'HTML')
        self.assertIn('summary', fastp.informs)

    def test_seqkit_fastp_se(self) -> None:
        """
        Tests fastp with SE input.
        """
        fastp = Fastp(self.camel)
        fastp.add_input_files({'FASTQ': [TestFastp.INPUT_FASTQ_PE[0]]})
        fastp.update_parameters(output_name='reads_illumina')
        fastp.run(self.running_dir)
        self.verify_output_files(fastp, 'FASTQ', nb_files=1)
        self.verify_output_files(fastp, 'JSON')
        self.verify_output_files(fastp, 'HTML')
        self.assertIn('summary', fastp.informs)


if __name__ == '__main__':
    Camel.get_instance()
    unittest.main()
