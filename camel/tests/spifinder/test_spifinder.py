import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.tools.pipelines.salmonella.spifinder import SPIFinder
from camel.app.tools.pipelines.salmonella.spifinderreporter import SPIFinderReporter


class TestSPIFinder(CamelTestSuite):
    """
    Tests the SPIFinder tool and its reporter.
    """
    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('salmonella')
    input_pe_reads = [test_file_dir / "SRR493330_1.fastq.gz", test_file_dir / "SRR493330_2.fastq.gz"]
    input_fasta_file = test_file_dir / 'assembly_filtered.fasta'
    output_file_spifinder_fastq = test_file_dir / 'output_spifinder_fastq.json'
    output_file_spifinder_fasta = test_file_dir / 'output_spifinder_fasta.json'
    db_path = Path('/db/pipelines/salmonella/spifinder/genomicepidemiology-spifinder_db-db102668b704')

    def test_spifinder_fastq(self) -> None:
        """
        Tests basic spifinder run on fastq input data.
        :return: None
        """
        spifinder = SPIFinder(self.camel)
        spifinder.add_input_files({
            'FASTQ_PE': [ToolIOFile(x) for x in self.input_pe_reads],
            'DIR': [ToolIODirectory(self.db_path)]
        })
        spifinder.run(self.running_dir)
        self.verify_output_files(spifinder, 'JSON')

    def test_spifinder_fasta(self) -> None:
        """
        Tests basic spifinder run on fasta input data.
        :return: None
        """
        spifinder = SPIFinder(self.camel)
        spifinder.add_input_files({
            'FASTA': [ToolIOFile(Path(TestSPIFinder.input_fasta_file))],
            'DIR': [ToolIODirectory(self.db_path)]
        })
        spifinder.run(self.running_dir)
        self.verify_output_files(spifinder, 'JSON')

    def test_spifinder_reporter(self) -> None:
        """
        Tests the spifinder reporter tool.
        :return: None
        """
        dummy_informs_spifinder_fasta = {'_name': 'test', }
        SnakemakeUtils.dump_object(dummy_informs_spifinder_fasta, Path('./informs_spifinder_fasta.dummy'))
        dummy_informs_spifinder_fastq = {'_name': 'test', }
        SnakemakeUtils.dump_object(dummy_informs_spifinder_fastq, Path('./informs_spifinder_fastq.dummy'))
        dummy_tsv_spifinder = Path('./dummy.tsv')
        dummy_tsv_spifinder.touch()
        dummy_tsv_spifinder_documentation = Path('./dummy2.tsv')
        dummy_tsv_spifinder_documentation.touch()
        spifinderreporter = SPIFinderReporter(self.camel)
        spifinderreporter.add_input_files({'TSV_output': [ToolIOFile(dummy_tsv_spifinder)],
                                           'TSV_documentation': [ToolIOFile(dummy_tsv_spifinder_documentation)],
                                           'JSON_FASTQ': [ToolIOFile(self.output_file_spifinder_fastq)],
                                           'JSON_FASTA': [ToolIOFile(self.output_file_spifinder_fasta)]})
        spifinderreporter.add_input_informs({'spifinder_fastq': dummy_informs_spifinder_fastq,
                                             'spifinder_fasta': dummy_informs_spifinder_fasta})
        spifinderreporter.run(self.running_dir)
        self.assertGreater(len(spifinderreporter.tool_outputs['VAL_HTML'][0].value.to_html()), 0)


if __name__ == '__main__':
    unittest.main()
