import filecmp
import unittest

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.muscle.muscle import Muscle


class TestMuscle(CamelTestSuite):
    """
    Tests the Muscle tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('muscle')

    def test_muscle(self) -> None:
        """
        Tests basic Muscle run.
        :return: None
        """
        muscle = Muscle(self.camel)
        muscle.add_input_files({
            'FASTA': [ToolIOFile(TestMuscle.test_file_dir / "unaligned.fasta")]
        })
        muscle.run(self.running_dir)
        self.verify_output_files(muscle, 'FASTA')

        compfile = TestMuscle.test_file_dir / 'aligned.fasta'
        self.assertTrue(filecmp.cmp(muscle.tool_outputs['FASTA'][0].path, compfile, shallow=False))

    def test_muscle_html(self) -> None:
        """
        Tests basic Muscle run with html output.
        :return: None
        """
        muscle = Muscle(self.camel)
        muscle.add_input_files({
            'FASTA': [ToolIOFile(TestMuscle.test_file_dir / "unaligned.fasta")]
        })
        muscle.update_parameters(html=None)
        muscle.run(self.running_dir)
        self.verify_output_files(muscle, 'HTML')

        compfile = TestMuscle.test_file_dir / 'aligned.html'
        self.assertTrue(filecmp.cmp(muscle.tool_outputs['HTML'][0].path, compfile, shallow=False))

    def test_muscle_clw(self) -> None:
        """
        Tests basic Muscle run with clw output.
        :return: None
        """
        muscle = Muscle(self.camel)
        muscle.add_input_files({
            'FASTA': [ToolIOFile(TestMuscle.test_file_dir / "unaligned.fasta")]
        })
        muscle.update_parameters(clw=None)
        muscle.run(self.running_dir)
        self.verify_output_files(muscle, 'CLW')

        compfile = TestMuscle.test_file_dir / 'aligned.clw'
        self.assertTrue(filecmp.cmp(muscle.tool_outputs['CLW'][0].path, compfile, shallow=False))

    def test_muscle_clw_strict(self) -> None:
        """
        Tests basic Muscle run with clw strict output.
        :return: None
        """
        muscle = Muscle(self.camel)
        muscle.add_input_files({
            'FASTA': [ToolIOFile(TestMuscle.test_file_dir / "unaligned.fasta")]
        })
        muscle.update_parameters(clwstrict=None)
        muscle.run(self.running_dir)
        self.verify_output_files(muscle, 'CLW')

        compfile = TestMuscle.test_file_dir / 'aligned.clws'
        self.assertTrue(filecmp.cmp(muscle.tool_outputs['CLW'][0].path, compfile, shallow=False))

    def test_muscle_msf(self) -> None:
        """
        Tests basic Muscle run with msf output.
        :return: None
        """
        muscle = Muscle(self.camel)
        muscle.add_input_files({
            'FASTA': [ToolIOFile(TestMuscle.test_file_dir / "unaligned.fasta")]
        })
        muscle.update_parameters(msf=None)
        muscle.run(self.running_dir)
        self.verify_output_files(muscle, 'MSF')

        compfile = TestMuscle.test_file_dir / 'aligned.msf'
        self.assertTrue(filecmp.cmp(muscle.tool_outputs['MSF'][0].path, compfile, shallow=False))

    def test_muscle_error(self) -> None:
        """
        Tests Muscle with broken fasta file.
        :return: None
        """
        muscle = Muscle(self.camel)
        muscle.add_input_files({
            'FASTA': [ToolIOFile(TestMuscle.test_file_dir / "bad.fasta")]
        })

        with self.assertRaises(ToolExecutionError):
            muscle.run(self.running_dir)


if __name__ == '__main__':
    unittest.main()
