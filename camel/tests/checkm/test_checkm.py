from pathlib import Path

from camel.app.camel import Camel
from camel.app.command.command import Command
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.checkm.checkm import CheckM
from camel.scripts.checkm.maincheckm import MainCheckM
from camel.tests import longRunningTest, resourceIntensiveTest


class TestCheckM(CamelTestSuite):
    """
    Tests the CheckM tool.
    """

    # Input files
    test_file_dir = CamelTestSuite.get_test_file_dir('checkm')
    input_fasta = test_file_dir / 'contigs_neisseria.fasta'

    def test_dependencies(self) -> None:
        """
        Tests if the tool dependencies are available.
        :return: None
        """
        checkm = CheckM(Camel.get_instance())
        for dependency in checkm.dependencies:
            command = Command(f'module load {dependency};')
            command.run_command(self.running_dir)
            self.assertEqual(command.returncode, 0, f"Dependency '{dependency}' cannot be loaded")

    @longRunningTest()
    @resourceIntensiveTest(reason='RAM usage')
    def test_checkm(self) -> None:
        """
        Tests the CheckM tool.
        :return: None
        """
        checkm = CheckM(Camel.get_instance())
        checkm.add_input_files({'FASTA': [ToolIOFile(TestCheckM.input_fasta)]})
        checkm.run(self.running_dir)
        self.assertIn('TSV', checkm.tool_outputs)
        self.assertGreater(Path(checkm.tool_outputs['TSV'][0].path).stat().st_size, 0)
        self.assertIn('results', checkm.informs)

    @longRunningTest()
    @resourceIntensiveTest(reason='RAM usage')
    def test_checkm_main_script(self) -> None:
        """
        Tests the CheckM main script.
        :return: None
        """
        path_report_out = self.running_dir / 'report' / 'report.html'
        checkv_main = MainCheckM([
            '--fasta', str(TestCheckM.input_fasta), TestCheckM.input_fasta.name,
            '--working-dir', str(self.running_dir),
            '--output-html', str(path_report_out),
            '--output-dir', str(path_report_out.parent)
        ])
        checkv_main.run()
        self.assertGreater(path_report_out.stat().st_size, 0)
