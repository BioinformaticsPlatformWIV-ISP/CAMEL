from camel.app.camel import Camel
from camel.app.command.command import Command
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.checkv.checkv import CheckV
from camel.scripts.checkv.maincheckv import MainCheckV
from camel.tests import resourceIntensiveTest, longRunningTest


class TestCheckV(CamelTestSuite):
    """
    Tests the CheckV tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('checkv')
    input_fasta = test_file_dir / 'contigs_hev.fasta'

    def test_dependencies(self) -> None:
        """
        Tests if the tool dependencies are available.
        :return: None
        """
        checkm = CheckV(Camel.get_instance())
        for dependency in checkm.dependencies:
            command = Command(f'module load {dependency};')
            command.run(self.running_dir)
            self.assertEqual(command.returncode, 0, f"Dependency '{dependency}' cannot be loaded")

    @longRunningTest()
    @resourceIntensiveTest(reason='RAM usage')
    def test_checkv_tool(self) -> None:
        """
        Tests the CheckV tool.
        :return: None
        """
        checkv = CheckV(Camel.get_instance())
        checkv.add_input_files({'FASTA': [ToolIOFile(TestCheckV.input_fasta)]})
        checkv.run(self.running_dir)
        self.assertIn('TSV_quality_summary', checkv.tool_outputs)

    @longRunningTest()
    @resourceIntensiveTest(reason='RAM usage')
    def test_checkv_main_script(self) -> None:
        """
        Tests the CheckV main script.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        checkv_main = MainCheckV([
            '--fasta', str(TestCheckV.input_fasta),
            '--working-dir', str(self.running_dir),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent)
        ])
        checkv_main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)
