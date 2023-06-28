import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.tools.integronfinder.integronfinder import IntegronFinder
from camel.app.tools.integronfinder.integronfinderreporter import IntegronFinderReporter
from camel.scripts.integronfinder.mainintegronfinder import MainIntegronFinder


class TestIntegronFinder(CamelTestSuite):
    """
    Tests for the IntegronFinder tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('integron_finder')
    FILE_FASTA = ToolIOFile(test_file_dir / 'GCF_002202175.1.fasta')
    FILE_FASTA_E_COLI = ToolIOFile(test_file_dir / 'NC_002695.2.fasta')

    def test_integron_finder(self) -> None:
        """
        Tests the IntegronFinder tool.
        :return: None
        """
        integron_finder = IntegronFinder(self.camel)
        integron_finder.add_input_files({'FASTA': [TestIntegronFinder.FILE_FASTA]})
        integron_finder.run(self.running_dir)
        self.verify_output_files(integron_finder, 'TSV')
        self.assertGreater(integron_finder.informs['nb_detected'], 0)

    def test_integron_finder_with_options(self) -> None:
        """
        Tests the IntegronFinder tool.
        :return: None
        """
        integron_finder = IntegronFinder(self.camel)
        integron_finder.add_input_files({'FASTA': [TestIntegronFinder.FILE_FASTA]})
        integron_finder.update_parameters(local_max=True, threads=8)
        integron_finder.run(self.running_dir)
        self.verify_output_files(integron_finder, 'TSV')
        self.assertGreater(integron_finder.informs['nb_detected'], 0)

    def test_integron_finder_reporter(self) -> None:
        """
        Tests the IntegronFinder tool with reporter.
        :return: None
        """
        # Run tool
        integron_finder = IntegronFinder(self.camel)
        integron_finder.add_input_files({'FASTA': [TestIntegronFinder.FILE_FASTA]})
        integron_finder.run(self.running_dir)
        self.verify_output_files(integron_finder, 'TSV')
        self.assertGreater(integron_finder.informs['nb_detected'], 0)

        # Run reporter
        reporter = IntegronFinderReporter(self.camel)
        reporter.add_input_files({'TSV': integron_finder.tool_outputs['TSV']})
        reporter.add_input_informs({'integron_finder': integron_finder.informs})
        reporter.run(self.running_dir)
        self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)

    def test_integron_finder_main_script(self) -> None:
        """
        Tests the IntegronFinder main script.
        """
        path_out_html = self.running_dir / 'out' / 'report.html'
        path_out_tsv = self.running_dir / 'out' / 'integrons.tsv'
        main_integron_finder = MainIntegronFinder([
            '--fasta', str(TestIntegronFinder.FILE_FASTA),
            '--working-dir', str(self.running_dir),
            '--output-tsv', str(path_out_tsv),
            '--output-html', str(path_out_html),
            '--output-dir', str(path_out_html.parent),
            '--threads', '4'
        ])
        main_integron_finder.run()

    def test_integron_finder_main_script_ecoli(self) -> None:
        """
        Tests the IntegronFinder main script on E. coli.
        """
        path_out_html = self.running_dir / 'out' / 'report.html'
        path_out_tsv = self.running_dir / 'out' / 'integrons.tsv'
        main_integron_finder = MainIntegronFinder([
            '--fasta', str(TestIntegronFinder.FILE_FASTA_E_COLI),
            '--working-dir', str(self.running_dir),
            '--output-tsv', str(path_out_tsv),
            '--output-html', str(path_out_html),
            '--output-dir', str(path_out_html.parent),
            '--threads', '4'
        ])
        main_integron_finder.run()


if __name__ == '__main__':
    unittest.main()
