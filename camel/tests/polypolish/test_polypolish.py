import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.polypolish.polypolish import Polypolish


class TestPolypolish(CamelTestSuite):
    """
    Initializes the Polypolish testing tool.
    """
    # Get test file and reference file directories
    test_file_dir = CamelTestSuite.get_test_file_dir('polypolish')

    # Create ToolIOFile input files
    FILE_SAM = ToolIOFile(test_file_dir / 'S1_ilmn.sam')
    FILE_FASTA_REF = ToolIOFile(test_file_dir / 'assembly.fasta')

    def test_polypolish(self) -> None:
        """
        Tests Polypolish.
        :return: None
        """
        polypolish = Polypolish(self.camel)
        polypolish.add_input_files({
            'SAM': [TestPolypolish.FILE_SAM],
            'FASTA': [TestPolypolish.FILE_FASTA_REF]})
        polypolish.run(self.running_dir)
        self.verify_output_files(polypolish, 'FASTA')


if __name__ == '__main__':
    unittest.main()
