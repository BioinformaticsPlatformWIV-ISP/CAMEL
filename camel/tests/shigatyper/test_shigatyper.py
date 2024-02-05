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
    Tests the ShigaTyper tool.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines')
    fastq_fwd = test_file_dir / 'Shigella-S17BD07654_1.fastq.gz'
    fastq_rev = test_file_dir / 'Shigella-S17BD07654_2.fastq.gz'

    sample_id = fastq_fwd.name.replace('_1.fastq.gz', '')

    def test_shigatyper(self) -> None:
        """
        Tests the ShigaTyper tool.
        :return: None
        """
        serotype_tool = ShigaTyper(self.camel)
        serotype_tool.add_input_files({'FASTQ_FWD': [ToolIOFile(Path(TestShigaTyper.fastq_fwd))],
                                       'FASTQ_REV': [ToolIOFile(Path(TestShigaTyper.fastq_rev))]})
        serotype_tool.run(self.running_dir)
        logging.info(f'Successfully processed: {TestShigaTyper.sample_id}')
        self.verify_output_files(serotype_tool, 'TSV')

    def test_shigatyper_reporter(self) -> None:
        """
        Tests the ShigaTyperReporter tool.
        :return: None
        """
        # Run the tool
        serotype_tool = ShigaTyper(self.camel)
        serotype_tool.add_input_files({'FASTQ_FWD': [ToolIOFile(Path(TestShigaTyper.fastq_fwd))],
                                       'FASTQ_REV': [ToolIOFile(Path(TestShigaTyper.fastq_rev))]})
        serotype_tool.run(self.running_dir)

        # Run the reporter
        reporter = ShigaTyperReporter(self.camel)
        reporter.add_input_files({'TSV': serotype_tool.tool_outputs['TSV'],
                                  'TSV_HITS': serotype_tool.tool_outputs['TSV_HITS']})
        reporter.add_input_informs({'shigatyper': serotype_tool.informs})
        reporter.run(self.running_dir)
        self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)

        # Save the report in a pickle
        path_html_io = self.running_dir / 'html.io'
        SnakemakeUtils.dump_object(reporter.tool_outputs['HTML'], path_html_io)
        logging.info(f'Report pickle saved to: {path_html_io}')


if __name__ == '__main__':
    unittest.main()
