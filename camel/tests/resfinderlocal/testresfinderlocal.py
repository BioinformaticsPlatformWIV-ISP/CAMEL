import argparse
import unittest
from pathlib import Path

import os
import tempfile

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.scripts.resfinderlocal.mainresfinderlocal import MainResFinderLocal


class TestResFinderLocal(unittest.TestCase):
    """
    Tests the ResFinder local tool.
    """

    camel = Camel()
    running_dir = None

    # Input files
    test_file_dir = Path(camel.config['testing']['testfiles_dir'])
    input_fasta = ToolIOFile(str(test_file_dir / 'workflows' / 'NC_002695.1.fasta'))
    input_fasta_galaxy = ToolIOFile(str(test_file_dir / 'workflows' / 'dataset_12.dat'))
    input_reads_no_hit = [ToolIOFile(str(test_file_dir / 'workflows' / 'ecoli_1.fastq')),
                          ToolIOFile(str(test_file_dir / 'workflows' / 'ecoli_2.fastq'))]
    input_reads_raw = [ToolIOFile(str(test_file_dir / 'gene_detection' / 'reads-ds_1P.fastq')),
                       ToolIOFile(str(test_file_dir / 'gene_detection' / 'reads-ds_2P.fastq'))]
    input_reads_raw_galaxy = [ToolIOFile(str(test_file_dir / 'workflows' / 'dataset_fwd_11.dat')),
                              ToolIOFile(str(test_file_dir / 'workflows' / 'dataset_rev_10.dat'))]
    input_db = '/db/gene_detection/ResFinder'

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(None, 'camel_', TestResFinderLocal.camel.config['temp_dir'])

    def __get_basic_arguments(self, report_path: str, detection_method: str) -> argparse.Namespace:
        """
        Returns the basic arguments for the main script.
        :param report_path: Report path
        :param detection_method: Detection method
        :return: Arguments
        """
        return argparse.Namespace(
            sample_name=None,
            fasta=None,
            fasta_name=None,
            fastq_pe=None,
            fastq_pe_names=None,
            resfinder_db=TestResFinderLocal.input_db,
            output_html=report_path,
            output_dir=os.path.dirname(report_path),
            trim_reads=True,
            min_percent_identity=90,
            min_percent_coverage=70,
            working_dir=self.running_dir,
            detection_method=detection_method,
            threads=4,
            report_include_fastq=False,
            kmers=55
        )

    def test_resfinder_local(self) -> None:
        """
        Tests the ResFinder local main script.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = self.__get_basic_arguments(output_file_report, 'blast')
        args.fasta = TestResFinderLocal.input_fasta.path
        main = MainResFinderLocal(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)

    def test_resfinder_local_spaces(self) -> None:
        """
        Tests the ResFinder local main script on an input file with spaces.
        :return: None
        """
        output_file_report = os.path.join(self.running_dir, 'report', 'report.html')
        args = self.__get_basic_arguments(output_file_report, 'blast')
        args.fasta = TestResFinderLocal.input_fasta.path
        args.fasta_name = 'reference with spaces.fasta'
        main = MainResFinderLocal(args)
        main.run()
        self.assertGreater(os.path.getsize(output_file_report), 0)


if __name__ == '__main__':
    unittest.main()
