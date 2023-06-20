import logging
import unittest
from pathlib import Path

from camel.app.components.html.htmlreport import HtmlReport
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.pipelines.neisseria.characterizeneisseriacapsulereporter import CharacterizeNeisseriaCapsuleReporter
from camel.app.tools.pipelines.neisseria.characterizeneisseriacapsule import CharacterizeNeisseriaCapsule
from camel.resources import CSS_STYLE


class TestCharacterizeNeisseriaCapsule(CamelTestSuite):
    """
    Tests the CharacterizeNeisseriaCapsule tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('characterizeneisseriacapsule')

    def test_characterizeneisseriacapsule(self) -> None:
        """
        Tests the CharacterizeNeisseriaCapsule tool.
        :return: None
        """
        # Run CharacterizeNeisseriaCapsule
        serogroup_tool = CharacterizeNeisseriaCapsule(self.camel)
        serogroup_tool.add_input_files({
                    'FASTA_dir': [ToolIODirectory(TestCharacterizeNeisseriaCapsule.test_file_dir)]})
        serogroup_tool.update_parameters(output_directory=f'{self.running_dir}')
        serogroup_tool.run(self.running_dir)
        ## TODO : method stops currently here due to warnings -> solve this
        logging.info(f'Successfully processed: {TestCharacterizeNeisseriaCapsule.test_file_dir}')

    def test_characterizeneisseriacapsulereporter(self) -> None:
        """
        Tests the CharacterizeNeisseriaCapsuleReporter tool.
        :return: None
        """
        # Run CharacterizeNeisseriaCapsule
        serogroup_tool = CharacterizeNeisseriaCapsule(self.camel)
        serogroup_tool.add_input_files({
                    'FASTA_dir': [ToolIODirectory(TestCharacterizeNeisseriaCapsule.test_file_dir)]})
        serogroup_tool.update_parameters(output_directory=f'{self.running_dir}')
        serogroup_tool.run(self.running_dir)
        logging.info(f'Successfully processed: {TestCharacterizeNeisseriaCapsule.test_file_dir}')

        # Run the reporter
        logging.info('Moved to reporter')
        reporter = CharacterizeNeisseriaCapsuleReporter(self.camel)
        ## TODO : update input files
        #reporter.add_input_files({'TSV': characterizeneisseriacapsule.tool_outputs['TSV']})
        reporter.add_input_files({'TSV': [ToolIOFile(Path('/scratch/temp/camel_12o24_we/serogroup/serogroup_predictions_1687254004.4534914.tab'))]})
        reporter.run(self.running_dir)

        # Save the report
        ## TODO : update HtmlReport input
        report = HtmlReport(Path('/scratch/temp/camel_12o24_we/html.report'), Path('/scratch/temp/camel_12o24_we/html.report').parent)
        report.initialize('Characterize Neisseria Capsule report on Sample', CSS_STYLE)
        report.add_html_object(reporter.tool_outputs['HTML'][0].value)
        report.save()

if __name__ == '__main__':
    unittest.main()
