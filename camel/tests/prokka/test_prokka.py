from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.fastx.fastqqualityfilter import FastqQualityFilter
from camel.app.tools.fastx.fastqqualitytrimmer import FastqQualityTrimmer
from camel.app.tools.prokka.prokka import Prokka


class TestProkka(CamelTestSuite):
    """
    Tests the Prokka tool.
    """

    # Input files
    input_fasta = ToolIOFile(CamelTestSuite.get_test_file_dir('prokka') / 'plasmid.fasta')

    def test_prokka(self) -> None:
        """
        Tests the Prokka tool.
        :return: None
        """
        prokka = Prokka(TestProkka.camel)
        prokka.add_input_files({'FASTA': [TestProkka.input_fasta]})
        prokka.run(self.running_dir)
        self.assertGreater(Path(prokka.tool_outputs['GFF'][0].path).stat().st_size, 0)
