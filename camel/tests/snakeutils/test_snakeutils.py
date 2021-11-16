import logging
import unittest
from pathlib import Path

import pkg_resources
import re

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.scripts.spatyping.mainspatyping import MainSpaTyping


class TestSpaTyping(CamelTestSuite):
    """
    Tests the spa typing tool.
    """
    # Input files
    test_snakefile = Path(pkg_resources.resource_filename('camel', 'tests/snakeutils/workflow_test.smk'))

    def test_run_snakefile(self) -> None:
        """
        Tests the run snakefile method.
        :return: None
        """
        config_data = {
            'working_dir': str(self.running_dir)
        }
        config_path = SnakePipelineUtils.generate_config_file(config_data, self.running_dir)
        SnakePipelineUtils.run_snakemake(str(TestSpaTyping.test_snakefile), config_path, [], self.running_dir)

    def test_run_snakefile_resources(self) -> None:
        """
        Tests the run snakefile method with the resources parameter set.
        :return: None
        """
        config_data = {
            'working_dir': str(self.running_dir)
        }
        config_path = SnakePipelineUtils.generate_config_file(config_data, self.running_dir)
        command = SnakePipelineUtils.run_snakemake(
            str(TestSpaTyping.test_snakefile), config_path, [], self.running_dir, resources={'RAM': 2, 'GPU': 4})
        self.assertIn('--resources', command.command, 'Resources parameter not added')
        self.assertIn('RAM', command.command, 'Resource not added')
        self.assertIn('GPU', command.command, 'Resource not added')


if __name__ == '__main__':
    unittest.main()
