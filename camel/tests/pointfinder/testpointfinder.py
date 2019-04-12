import argparse
import unittest

import os
import tempfile

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.scripts.pointfinder.mainpointfinder import MainPointFinder


class TestPointFinder(unittest.TestCase):
    """
    Tests the PointFinder tool.
    """
    camel = Camel()
    running_dir = None

    # Input files
    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'])
    input_fasta_file = ToolIOFile(os.path.join(test_file_dir, 'pointfinder', 'ref_ecoli.fasta'))
    input_fasta_file_salm = ToolIOFile(os.path.join(test_file_dir, 'pointfinder', 'salmonella_lt2_ref.fasta'))
    input_fastq_raw_galaxy = [ToolIOFile('/data/camel/testdata/workflows/dataset_fwd_11.dat'),
                              ToolIOFile('/data/camel/testdata/workflows/dataset_rev_10.dat')]

    def setUp(self) -> None:
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(prefix='camel_', dir=TestPointFinder.camel.config['temp_dir'])

    def test_pointfinder(self) -> None:
        """
        Tests the PointFinder main script
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            sample_name='test_sample',
            fasta=self.input_fasta_file.path,
            fasta_name=os.path.basename(self.input_fasta_file.basename),
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            working_dir=self.running_dir,
            species='e.coli'
        )
        main = MainPointFinder(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_pointfinder_with_pubmed_link(self) -> None:
        """
        Tests the PointFinder main script with a link to PubMed.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            sample_name='test_sample',
            fasta=self.input_fasta_file_salm.path,
            fasta_name=os.path.basename(self.input_fasta_file_salm.basename),
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            working_dir=self.running_dir,
            species='e.coli'
        )
        main = MainPointFinder(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_pointfinder_fastq_input(self) -> None:
        """
        Tests the PointFinder main script with FASTQ input.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = argparse.Namespace(
            sample_name='test_sample',
            fasta=None,
            fasta_name=None,
            fastq_pe=[f.path for f in self.input_fastq_raw_galaxy],
            fastq_pe_names=['my-sample_R1.fastq', 'my-sample_R2.fastq'],
            output_html=output_file_report,
            output_dir=os.path.dirname(output_file_report),
            working_dir=self.running_dir,
            species='e.coli',
            trim_reads=False,
            kmers=None,
            threads=4
        )
        main = MainPointFinder(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)


if __name__ == '__main__':
    unittest.main()
