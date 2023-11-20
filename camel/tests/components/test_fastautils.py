import unittest

from camel.app.components.files.fastautils import FastaUtils
from camel.app.components.testing.cameltestsuite import CamelTestSuite


class TestFastaUtils(CamelTestSuite):
    """
    Tests the Fasta utils module.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('components')


    def test_convert_fasta_to_fastq(self) -> None:
        """
        Test the function to convert a FASTA file to a FASTQ file
        :return: None
        """
        input_file = TestFastaUtils.test_file_dir / 'toy.fasta'
        output_file = self.running_dir / f"{input_file.stem}.fastq"
        FastaUtils.convert_fasta_to_fastq(input_file, output_file)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
