import unittest
from pathlib import Path
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.tools.pipelines.salmonella.sistr import Sistr
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile


class TestSistr(CamelTestSuite):
    """
    Tests the Sistr tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('salmonella')
    input_fasta_file = test_file_dir / 'assembly_filtered.fasta'

    def test_sistr(self) -> None:
        """
        Tests basic Sistr run.
        :return: None
        """
        sistrtool = Sistr(self.camel)
        sistrtool.add_input_files({
            'FASTA': [ToolIOFile(Path(TestSistr.input_fasta_file))],
            'DIR': [ToolIODirectory(Path('/db/SISTR/1.1.1/data'))]
        })
        sistrtool.run(self.running_dir)
        self.verify_output_files(sistrtool, 'JSON')


if __name__ == '__main__':
    unittest.main()
