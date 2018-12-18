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
        self.assertEquals(FastqUtils.get_sample_name('/data/temp/Z4686_S31_L001_R1_001.fastq.gz'), 'Z4686')

    def test_get_sample_name_miseq_fmt_no_s(self) -> None:
        """
        Tests the get sample name function for MiSEQ format without the sample number.
        :return: None
        """
        self.assertEquals(FastqUtils.get_sample_name('/data/temp/Z4686_L001_R1_001.fastq.gz'), 'Z4686')

    def test_get_sample_name_miseq_fmt_no_s_no_l(self) -> None:
        """
        Tests the get sample name function for MiSEQ format without the sample number and lane number.
        :return: None
        """
        self.assertEquals(FastqUtils.get_sample_name('/data/temp/Z4686_R1_001.fastq.gz'), 'Z4686')

    def test_get_sample_name_simple_fmt(self) -> None:
        """
        Tests the get sample name function for the simple format.
        :return: None
        """
        self.assertEquals(FastqUtils.get_sample_name('/data/temp/Z4686_1.fastq.gz'), 'Z4686')

    def test_get_sample_name_invalid_fmt(self) -> None:
        """
        Tests the get sample name function for an invalid filename.
        :return: None
        """
        with self.assertRaises(ValueError):
            FastqUtils.get_sample_name('/data/temp/invalid_name.fastq.gz')


if __name__ == '__main__':
    unittest.main()
