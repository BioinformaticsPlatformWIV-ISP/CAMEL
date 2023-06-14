import logging
import unittest
import os
from pathlib import Path

from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.components.testing.cameltestsuite import CamelTestSuite
#from camel.app.tools.pipelines.neisseria.characterizeneisseriacapsulereporter import CharacterizeNeisseriaCapsuleReporter
from camel.app.tools.pipelines.neisseria.characterizeneisseriacapsule import CharacterizeNeisseriaCapsule


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
        my_tool = CharacterizeNeisseriaCapsule(self.camel)
        my_tool.add_input_files({
                    'FASTA_dir': [ToolIODirectory(TestCharacterizeNeisseriaCapsule.test_file_dir)]
                })
        my_tool.update_parameters(output_directory=Path(f'{self.running_dir}'))
        my_tool.run()
        logging.info(f'Successfully processed: {TestCharacterizeNeisseriaCapsule.test_file_dir}')

if __name__ == '__main__':
    unittest.main()
