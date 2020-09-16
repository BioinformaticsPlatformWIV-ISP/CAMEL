from pathlib import Path

from camel.app.camel import Camel
from camel.app.command.command import Command
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.checkm.checkm import CheckM
from camel.tests import longRunningTest


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
    def test_checkm(self) -> None:
        """
        Tests the checkm tool.
        :return: None
        """
        checkm = CheckM(Camel.get_instance())
        checkm.add_input_files({'FASTA': [ToolIOFile(TestCheckM.input_fasta)]})
        checkm.run(self.running_dir)
        self.assertIn('TSV', checkm.tool_outputs)
        self.assertGreater(Path(checkm.tool_outputs['TSV'][0].path).stat().st_size, 0)
        self.assertIn('results', checkm.informs)
