import unittest
from pathlib import Path

from camel.app.components.testing.cameltestsuite import CamelTestSuite
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.star.staralign import StarAlign
from camel.app.tools.star.starindex import StarIndex


class TestStar(CamelTestSuite):
    """
    Tests the STAR tool.
    """

    test_file_dir = CamelTestSuite.get_test_file_dir("star")
    FILE_REF_GENOME = ToolIOFile(test_file_dir / "reference.fasta")
    FILE_FASTQ = ToolIOFile(test_file_dir / "reads.fastq")
    INDEX_DIR = ToolIODirectory(test_file_dir / "index")

    def test_star_index(self) -> None:
        """
        Tests StarIndex.
        :return: None
        """
        # set index path explicitly to avoid Permission Denied error on /testdata/camel/star
        path_index_dir = Path(self.running_dir) / "STAR_index"
        path_index_dir.mkdir(parents=True, exist_ok=True)

        star_index = StarIndex(self.camel)
        star_index.add_input_files({
            'FASTA': [TestStar.FILE_REF_GENOME],
            'INDEX_DIR': [ToolIODirectory(path_index_dir)]
        })
        star_index.run(self.running_dir)
        self.assertTrue('INDEX_DIR' in star_index.tool_outputs, "No index generated")
        output_dir = Path(str(star_index.tool_outputs['INDEX_DIR'][0]))
        self.assertTrue(any(output_dir.iterdir()))

    def test_star_align(self) -> None:
        """
        Tests StarAlign (single ends).
        :return: None
        """
        star_align = StarAlign(self.camel)
        star_align.add_input_files({
            'FASTQ': [TestStar.FILE_FASTQ],
            'INDEX_DIR': [TestStar.INDEX_DIR]
        })
        star_align.run(self.running_dir)
        self.assertTrue('ALIGNMENT' in star_align.tool_outputs, "No SAM output generated")
        output_file = Path(star_align.tool_outputs['ALIGNMENT'][0].path)
        self.assertTrue(output_file.exists())
        self.assertGreater(output_file.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
