import unittest

import os
import tempfile

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.fastx.fastqqualityfilter import FastqQualityFilter
from camel.app.tools.fastx.fastqqualitytrimmer import FastqQualityTrimmer


class TestFastx(unittest.TestCase):
    """
    Tests the fastx tools.
    """
    camel = Camel()
    running_dir = None

    # Input files
    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'])
    input_fastq_file = ToolIOFile(os.path.join(test_file_dir, 'fastx', 'ERR2019997-iontorrent-ds.fastq'))

    def setUp(self) -> None:
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(prefix='camel_', dir=TestFastx.camel.config['temp_dir'])

    def test_quality_filter(self) -> None:
        """
        Tests the fastx quality filter.
        :return: None
        """
        q_filter = FastqQualityFilter(TestFastx.camel)
        q_filter.update_parameters(min_quality=30)
        q_filter.add_input_files({'FASTQ': [TestFastx.input_fastq_file]})
        q_filter.run(self.running_dir)
        self.assertGreater(TestFastx.input_fastq_file.size, q_filter.tool_outputs['FASTQ'][0].size)
        self.assertGreater(q_filter.informs['input_reads'], q_filter.informs['output_reads'])

    def test_quality_trimmer(self) -> None:
        """
        Tests the fastx quality trimmer.
        :return: None
        """
        q_trimmer = FastqQualityTrimmer(TestFastx.camel)
        q_trimmer.add_input_files({'FASTQ': [TestFastx.input_fastq_file]})
        q_trimmer.run(self.running_dir)
        self.assertGreater(TestFastx.input_fastq_file.size, q_trimmer.tool_outputs['FASTQ'][0].size)
        self.assertGreater(q_trimmer.informs['input_reads'], q_trimmer.informs['output_reads'])
