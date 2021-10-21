import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.bwa.bwaindex import BWAIndex
from camel.app.tools.bwa.bwamap import BWAMap

class TestBWA(CamelTestSuite):
    """
    Tests the BWA tool suite.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('bwa')

    def test_bwa_index(self) -> None:
        """
        Test BWAIndex
        :return: None
        """
        bwa_index = BWAIndex(self.camel)
        bwa_index.add_input_files({
            'FASTA_REF': [ToolIOFile(TestBWA.test_file_dir / 'Homo_sapiens_assembly38_chr22.fasta')],
        })
        bwa_index.run(self.running_dir)
        self.assertTrue('INDEX_GENOME_PREFIX' in bwa_index.tool_outputs, "No INDEX_GENOME_PREFIX output generated")
        output_file = Path(bwa_index.tool_outputs['INDEX_GENOME_PREFIX'][0].value)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

    def test_bwa_map(self) -> None:
        """
        Test BWAMap
        :return: None
        """
        bwa_map = BWAMap(self.camel)
        bwa_map.add_input_files({
            'FASTQ_PE': [
                ToolIOFile(TestBWA.test_file_dir / 'NA12877_R1.fastq'),
                ToolIOFile(TestBWA.test_file_dir / 'NA12877_R2.fastq')
            ],
            'INDEX_GENOME_PREFIX': [ToolIOValue(TestBWA.test_file_dir / 'Homo_sapiens_assembly38_chr22.fasta')],
        })
        bwa_map.run(self.running_dir)
        self.assertTrue('SAM' in bwa_map.tool_outputs, "No SAM output generated")
        output_file = Path(bwa_map.tool_outputs['SAM'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)

if __name__ == '__main__':
    unittest.main()