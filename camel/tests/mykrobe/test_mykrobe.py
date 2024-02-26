import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.tools.pipelines.salmonella.genotyphi import Genotyphi
from camel.app.tools.pipelines.salmonella.genotyphireporter import GenotyphiReporter


class TestGenotyphi(CamelTestSuite):
    """
    Tests the genotyphi tool.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('salmonella')
    input_pe_reads = [test_file_dir / "SRR493330_1.fastq.gz", test_file_dir / "SRR493330_2.fastq.gz"]
    genotyphi_output_csv = test_file_dir / 'genotyphi.csv'

    def test_genotyphi(self) -> None:
        """
        Tests basic genotyphi run.
        :return: None
        """
        genotyphitool = Genotyphi(self.camel)
        genotyphitool.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in self.input_pe_reads],
            'DIR': [ToolIODirectory(Path('/db/pipelines/salmonella/mykrobe/20220331'))]
        })
        genotyphitool.run(self.running_dir)
        self.verify_output_files(genotyphitool, 'CSV')

    def test_genotyphi_reporter(self) -> None:
        """
        Tests genotyphi reporter tool
        :return: None
        """
        dummy_informs_genotyphi = {'_name': 'test', }
        SnakemakeUtils.dump_object(dummy_informs_genotyphi, Path('./informs.dummy'))
        dummy_tsv_genotyphi = Path('./dummy.tsv')
        dummy_tsv_genotyphi.touch()
        genotyphireporter = GenotyphiReporter(self.camel)
        genotyphireporter.add_input_files({'TSV_output': [ToolIOFile(dummy_tsv_genotyphi)],
                                           'CSV': [ToolIOFile(self.genotyphi_output_csv)]})
        genotyphireporter.add_input_informs({'genotyphi': dummy_informs_genotyphi})
        genotyphireporter.run(self.running_dir)
        self.assertGreater(len(genotyphireporter.tool_outputs['VAL_HTML'][0].value.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
