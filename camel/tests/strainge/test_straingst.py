import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.strainge.straingstrun import StrainGSTRun
from camel.app.tools.strainge.straingstkmerize import StrainGSTKmerize


class TestStrainGST(CamelTestSuite):
    """
    Initializes this testing tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('strainge')
    FILE_FASTQ = ToolIOFile(test_file_dir / 'bsubtilis.fastq.gz')
    FILE_HDF5 = ToolIOFile(test_file_dir / 'bsubtilis.hdf5')
    DB_HDF5 = ToolIOFile(test_file_dir / 'small_bacillus_db.hdf5')

    def test_straingst_kmerize(self) -> None:
        """
        Testing StrainGST kmerize on ont sequencing data.
        :return: None
        """
        straingst_kmerize = StrainGSTKmerize(self.camel)
        straingst_kmerize.add_input_files({'FASTQ': [TestStrainGST.FILE_FASTQ]})
        straingst_kmerize.run(self.running_dir)
        self.verify_output_files(straingst_kmerize, 'HDF5')

    def test_straingst_run(self) -> None:
        """
        Testing StrainGST run on ONT sequencing data and a small B. subtilis database.
        :return: None
        """
        straingst_run = StrainGSTRun(self.camel)
        straingst_run.add_input_files({'HDF5': [TestStrainGST.FILE_HDF5], 'DB_HDF5': [TestStrainGST.DB_HDF5]})
        straingst_run.run(self.running_dir)
        self.verify_output_files(straingst_run, 'TSV_STATS')


if __name__ == '__main__':
    unittest.main()
