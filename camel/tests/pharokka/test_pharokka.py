import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.pharokka.pharokka import Pharokka
from camel.app.tools.pharokka.pharokkareporter import PharokkaReporter
from camel.tests import minOSVersion


class TestPharokka(CamelTestSuite):
    """
    Tests the Pharokka tool.
    """
    test_file_dir = CamelTestSuite.get_test_file_dir('pharokka')
    # fasta = test_file_dir / 'S_20_721_prophage.fasta'
    fasta = test_file_dir / 'pBAD33.fasta'

    @minOSVersion('jammy')
    def test_pharokka(self) -> None:
        """
        Tests the Pharokka tool.
        :return: None
        """
        pharokka = Pharokka(self.camel)
        pharokka.add_input_files({'FASTA': [ToolIOFile(Path(TestPharokka.fasta))]})
        pharokka.run(self.running_dir)
        self.verify_output_files(pharokka, 'GBK')
        self.verify_output_files(pharokka, 'CARD')
        self.verify_output_files(pharokka, 'VFDB')

    @minOSVersion('jammy')
    def test_pharokka_reporter(self) -> None:
        """
        Tests the PharokkaReporter tool.
        :return: None
        """
        # Run the tool
        pharokka = Pharokka(self.camel)
        pharokka.add_input_files({'FASTA': [ToolIOFile(Path(TestPharokka.fasta))]})
        pharokka.run(self.running_dir)

        # Run the reporter
        reporter = PharokkaReporter(self.camel)
        reporter.add_input_files({'STATS': pharokka.tool_outputs['STATS'],
                                  'CARD': pharokka.tool_outputs['CARD'],
                                  'VFDB': pharokka.tool_outputs['VFDB'],
                                  'INPHARED': pharokka.tool_outputs['INPHARED'],
                                  'GBK': pharokka.tool_outputs['GBK']
                                  })
        reporter.add_input_informs({'pharokka': pharokka.informs})
        reporter.run(self.running_dir)
        self.assertGreater(len(reporter.tool_outputs['HTML'][0].value.to_html()), 0)
        CamelTestSuite.export_report_section(reporter.tool_outputs['HTML'][0].value, self.running_dir / 'report')


if __name__ == '__main__':
    unittest.main()
