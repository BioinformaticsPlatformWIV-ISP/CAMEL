from pathlib import Path

from camel.app.camel import Camel
from camel.app.command.command import Command
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.confindr.confindr import ConFindr
from camel.app.tools.confindr.confindrreporter import ConFindrReporter
from camel.scripts.confindr.mainconfindr import MainConFindr


class TestConFindr(CamelTestSuite):
    """
    Tests the ConFindr tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('pipelines')
    input_pe_reads = [
        test_file_dir / 'Neisseria-2011-006_S6-ds_1.fastq.gz',
        test_file_dir / 'Neisseria-2011-006_S6-ds_2.fastq.gz'
    ]

    def test_dependencies(self) -> None:
        """
        Tests if the tool dependencies are available.
        :return: None
        """
        confindr = ConFindr(Camel.get_instance())
        for dependency in confindr.dependencies:
            command = Command(f'module load {dependency};')
            command.run_command(self.running_dir)
            self.assertEqual(command.returncode, 0, f"Dependency '{dependency}' cannot be loaded")

    def test_confindr_pe(self) -> None:
        """
        Tests confindr with PE input.
        :return: None
        """
        confindr = ConFindr(Camel.get_instance())
        confindr.add_input_files({
            'FASTQ_GZ_PE': [ToolIOFile(x) for x in TestConFindr.input_pe_reads]
        })
        confindr.run(self.running_dir)
        self.assertIn('CSV', confindr.tool_outputs)
        self.assertGreater(Path(confindr.tool_outputs['CSV'][0].path).stat().st_size, 0)
        self.assertIn('ContamStatus', confindr.informs)

    def test_confindr_se(self) -> None:
        """
        Tests confindr with PE input.
        :return: None
        """
        confindr = ConFindr(Camel.get_instance())
        confindr.add_input_files({
            'FASTQ_GZ_SE': [ToolIOFile(TestConFindr.input_pe_reads[0])]
        })
        confindr.run(self.running_dir)
        self.assertIn('CSV', confindr.tool_outputs)
        self.assertGreater(Path(confindr.tool_outputs['CSV'][0].path).stat().st_size, 0)
        self.assertIn('ContamStatus', confindr.informs)

    def test_confindr_se_with_report(self) -> None:
        """
        Tests confindr with PE input.
        :return: None
        """
        confindr = ConFindr(Camel.get_instance())
        confindr.add_input_files({
            'FASTQ_GZ_SE': [ToolIOFile(TestConFindr.input_pe_reads[0])]
        })
        confindr.run(self.running_dir)
        self.assertIn('CSV', confindr.tool_outputs)
        self.assertGreater(Path(confindr.tool_outputs['CSV'][0].path).stat().st_size, 0)
        self.assertIn('ContamStatus', confindr.informs)

        reporter = ConFindrReporter(Camel.get_instance())
        reporter.add_input_informs({'confindr': confindr.informs})
        reporter.run(self.running_dir)
        self.assertIn('HTML', reporter.tool_outputs)
        self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)

    def test_confindr_main_script(self) -> None:
        """
        Tests the ConFinder main script.
        :return: None
        """
        dir_out = self.running_dir / 'out'
        dir_out.mkdir()
        confindr_main = MainConFindr([
            '--fastq-se', str(TestConFindr.input_pe_reads[0]),
            '--working-dir', str(self.running_dir),
            '--output-html', str(dir_out / 'report.html'),
            '--output-dir', str(dir_out),
            '--quality-cutoff', '15',
            '--base-cutoff', '5',
            '--base-percentage-cutoff', '10',
            '--min-matching-hashes', '200',
        ])
        confindr_main.run()
