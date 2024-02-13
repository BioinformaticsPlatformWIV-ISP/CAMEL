import logging
import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.tools.pipelines.shigella.shigatyper import ShigaTyper
from camel.app.tools.pipelines.shigella.shigatyperreporter import ShigaTyperReporter


class TestShigaTyper(CamelTestSuite):
    """
    Tests the ShigaTyper tool and the related reporter tool.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines')
    fastq_fwd = test_file_dir / 'Shigella-S17BD07654_1.fastq.gz'
    fastq_rev = test_file_dir / 'Shigella-S17BD07654_2.fastq.gz'

    def test_shigatyper(self) -> None:
        """
        Tests the ShigaTyper tool.
        :return: None
        """
        shigatyper = ShigaTyper(self.camel)
        shigatyper.add_input_files({'FASTQ_PE': [ToolIOFile(Path(TestShigaTyper.fastq_fwd)),
                                                 ToolIOFile(Path(TestShigaTyper.fastq_rev))]})
        shigatyper.run(self.running_dir)
        self.verify_output_files(shigatyper, 'TSV')
        self.verify_output_files(shigatyper, 'TSV_HITS')

    def test_shigatyper_reporter(self) -> None:
        """
        Tests the ShigaTyperReporter tool.
        :return: None
        """
        # Run the tool
        shigatyper = ShigaTyper(self.camel)
        shigatyper.add_input_files({'FASTQ_PE': [ToolIOFile(Path(TestShigaTyper.fastq_fwd)),
                                                 ToolIOFile(Path(TestShigaTyper.fastq_rev))]})
        shigatyper.run(self.running_dir)

        # Run the reporter
        reporter = ShigaTyperReporter(self.camel)
        reporter.add_input_files({'TSV': shigatyper.tool_outputs['TSV'],
                                  'TSV_HITS': shigatyper.tool_outputs['TSV_HITS']})
        reporter.add_input_informs({'shigatyper': shigatyper.informs})
        reporter.run(self.running_dir)
        self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)

        # Save the report in a pickle
        path_html_io = self.running_dir / 'html.io'
        SnakemakeUtils.dump_object(reporter.tool_outputs['HTML'], path_html_io)
        logging.info(f'Report pickle saved to: {path_html_io}')


if __name__ == '__main__':
    unittest.main()
