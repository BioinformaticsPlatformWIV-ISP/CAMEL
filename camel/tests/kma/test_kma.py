import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.bcftools.bcftoolscsq import BcftoolsCsq
from camel.app.tools.bcftools.bcftoolsfilter import BcftoolsFilter
from camel.app.tools.kma.kma import KMA


class TestKMA(CamelTestSuite):
    """
    Tests the KMA tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('gene_detection')
    input_reads_raw = [test_file_dir / 'reads-ds_1P.fastq', test_file_dir / 'reads-ds_2P.fastq']
    input_gene_detection_db = test_file_dir / 'db' / 'kma' / 'arg-annot-clustered_80'

    def test_kma_pe(self) -> None:
        """
        Tests KMA with paired-end input.
        :return: None
        """
        kma = KMA(self.camel)
        kma.add_input_files({
            'FASTQ_PE': [ToolIOFile(file_) for file_ in TestKMA.input_reads_raw],
            'DB': [ToolIOValue(str(TestKMA.input_gene_detection_db))]
        })
        kma.run(self.running_dir)
        self.assertTrue('TSV' in kma.tool_outputs, "No TSV output generated")
        output_file = Path(kma.tool_outputs['TSV'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_kma_se(self) -> None:
        """
        Tests KMA with single-end input.
        :return: None
        """
        kma = KMA(self.camel)
        kma.add_input_files({
            'FASTQ_SE': [ToolIOFile(TestKMA.input_reads_raw[0])],
            'DB': [ToolIOValue(str(TestKMA.input_gene_detection_db))]
        })
        kma.run(self.running_dir)
        self.assertTrue('TSV' in kma.tool_outputs, "No TSV output generated")
        output_file = Path(kma.tool_outputs['TSV'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)



if __name__ == '__main__':
    unittest.main()
