import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.medaka.medakaconsensus import MedakaConsensus
from camel.app.tools.medaka.medakastitch import MedakaStitch


class TestMedaka(CamelTestSuite):
    """
    Tests the Medaka tool suite.
    """
    # Get test file and reference file directories
    test_file_dir = CamelTestSuite.get_test_file_dir('medaka')

    # Create ToolIOFile input files
    FILE_BAM = ToolIOFile(test_file_dir / 'calls_to_draft_subsampled.bam')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'contig_1.fasta')
    FILE_HDF = ToolIOFile(test_file_dir / 'consensus_probs_subsample.hdf')

    def test_medaka_consensus(self) -> None:
        """
        Tests MedakaConsensus.
        :return: None
        """
        consensus = MedakaConsensus(self.camel)
        consensus.add_input_files({'BAM': [TestMedaka.FILE_BAM]})
        consensus.run(self.running_dir)
        self.verify_output_files(consensus, 'HDF')

    def test_medaka_stitch(self) -> None:
        """
        Tests MedakaStitch.
        :return: None
        """
        stitch = MedakaStitch(self.camel)
        stitch.add_input_files({
            'HDF': [TestMedaka.FILE_HDF],
            'FASTA': [TestMedaka.FILE_FASTA_REF]
        })
        stitch.run(self.running_dir)
        self.verify_output_files(stitch, 'FASTA')


if __name__ == '__main__':
    unittest.main()
