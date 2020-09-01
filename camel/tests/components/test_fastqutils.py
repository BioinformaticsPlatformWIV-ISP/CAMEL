import unittest

from camel.app.components.files.fastqutils import FastqUtils


class TestFastqUtils(unittest.TestCase):
    """
    Tests the Fastq utils module.
    """

    def test_get_sample_name_miseq_fmt(self) -> None:
        """
        Tests the get sample name function for MiSEQ format.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name('/data/temp/Z4686_S31_L001_R1_001.fastq.gz'), 'Z4686')

    def test_get_sample_name_miseq_fmt_no_s(self) -> None:
        """
        Tests the get sample name function for MiSEQ format without the sample number.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name('/data/temp/Z4686_L001_R1_001.fastq.gz'), 'Z4686')

    def test_get_sample_name_miseq_fmt_no_s_lowercase(self) -> None:
        """
        Tests the get sample name function for MiSEQ format without the sample number.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name('/data/temp/Z4686_L001_R1_001.fastq.gz'), 'Z4686')

    def test_get_sample_name_miseq_fmt_no_s_no_l(self) -> None:
        """
        Tests the get sample name function for MiSEQ format without the sample number and lane number.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name('/data/temp/Z4686_R1_001.fastq.gz'), 'Z4686')

    def test_get_sample_name_simple_fmt(self) -> None:
        """
        Tests the get sample name function for the simple format.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name('/data/temp/Z4686_1.fastq.gz'), 'Z4686')

    def test_get_sample_name_simple_fmt_alt(self) -> None:
        """
        Tests the get sample name function for the simple format.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name('/data/temp/S18BD02705_R1.fastq.gz'), 'S18BD02705')

    def test_get_sample_name_invalid_fmt(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        with self.assertRaises(ValueError):
            FastqUtils.get_sample_name('/data/temp/invalid_name.fastq.gz')

    def test_get_sample_name_with_dots(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name('/data/temp/my-sample.1.fastq.gz'), 'my-sample')

    def test_get_sample_name_fq_ext(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name('/data/temp/my-sample_1.fq'), 'my-sample')

    def test_get_sample_name_fq_ext_gzipped(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name('/data/temp/my-sample_1.fq.gz'), 'my-sample')

    def test_get_sample_name_with_p(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name('/data/temp/my-sample_1P.fastq'), 'my-sample')

    def test_get_sample_name_with_underscore(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name('/data/temp/my_sample_1.fastq'), 'my_sample')

    def test_get_sample_name_se(self) -> None:
        """
        Tests the get sample name function for a single end sample name.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(
            '/data/temp/my_sample.fastq', FastqUtils.PATTERN_FQ_SE), 'my_sample')

    def test_get_sample_name_se_gzipped(self) -> None:
        """
        Tests the get sample name function for a single end sample name.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(
            '/data/temp/my_sample.fastq.gz', FastqUtils.PATTERN_FQ_SE), 'my_sample')

    def test_get_sample_name_parentheses(self) -> None:
        """
        Tests the get sample name function that contains parentheses.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name('/data/temp/UI-586(SRR7648453)_1.fastq.gz'), 'UI-586SRR7648453')

    def test_get_sample_name_miseq(self) -> None:
        """
        Tests the get sample name function that contains parentheses.
        :return: None
        """
        self.assertEqual(FastqUtils.get_sample_name(
            '/data/temp/S20BD03018_S16_L001_R1_001.fastq.gz'), 'S20BD03018')


if __name__ == '__main__':
    unittest.main()
