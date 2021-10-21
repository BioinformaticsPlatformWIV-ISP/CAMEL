import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.kraken.kraken2 import Kraken2


class TestKraken2(CamelTestSuite):
    """
    Tests the Kraken2 tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('kraken')
    input_pe_reads = [test_file_dir / f'lm1_1.fastq', test_file_dir / f'lm1_2.fastq']

    def test_kraken2_paired(self) -> None:
        """
        Tests krqken2 with paired-end input.
        :return: None
        """
        kraken2 = Kraken2(TestKraken2.camel)
        kraken2.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in self.input_pe_reads],
            'DB': [ToolIODirectory(Path('/db/kraken2/latest/abfhpv/'))]
        })
        kraken2.run(self.running_dir)
        self.assertGreater(Path(kraken2.tool_outputs['TSV'][0].path).stat().st_size, 0)
        self.assertGreater(Path(kraken2.tool_outputs['TSV_report'][0].path).stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
