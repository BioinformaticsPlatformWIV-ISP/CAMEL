import logging
import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.pipelines.neisseria.gmats import GMats


class TestGMats(CamelTestSuite):
    """
    Tests the gMATS tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir('gmats')
    DB = Path(CamelTestSuite.camel.config['db_root'], 'pipelines', 'neisseria', 'gMATS_DB.txt')

    def test_gmats(self) -> None:
        """
        Tests the gMATS tool.
        :return: None
        """
        for tsv_path in TestGMats.test_file_dir.glob('*.tsv'):
            # Create the output directory
            output_dir = self.running_dir / f'{tsv_path.stem}'
            output_dir.mkdir(parents=True)
            logging.info(f'Saving output in: {output_dir}')

            # Run gMATS
            gmats = GMats(self.camel)
            gmats.add_input_files({
                'TSV': [ToolIOFile(tsv_path)],
                'DB': [ToolIOFile(TestGMats.DB)]
            })
            gmats.update_parameters(output_directory=str(output_dir))
            gmats.run(self.running_dir)
            self.verify_output_files(gmats, 'TSV')
            logging.info(f'Successfully processed: {tsv_path}')


if __name__ == '__main__':
    unittest.main()
