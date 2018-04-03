import unittest

import logging
import os
import tempfile
import yaml

from app.camel import Camel
from app.io.tooliofile import ToolIOFile
from app.snakemake.snakemakeutils import SnakemakeUtils
from resources.snakefile import SNAKEFILE_GENE_DETECTION
from resources.snakefile.gene_detection import OUTPUT_GENE_DETECTION_REPORT, INPUT_GENE_DETECTION_FASTA, \
    OUTPUT_GENE_DETECTION_REPORT_EMPTY, INPUT_GENE_DETECTION_FASTQ_PE
from run_gene_detection import CONFIG_TEMPLATE
from app.snakemake.snakemaketestutils import SnakemakeTestUtils


class TestGeneDetection(unittest.TestCase):
    """
    Tests the gene detection workflow.
    """

    camel = Camel()
    running_dir = None

    test_file_dir = os.path.join(camel.config['testing']['testfiles_dir'], 'gene_detection')

    FILE_FASTA = ToolIOFile(os.path.join(test_file_dir, 'contigs.fasta'))
    FILES_FASTQ_PE = [ToolIOFile(os.path.join(test_file_dir, 'reads-ds_1P.fastq')),
                      ToolIOFile(os.path.join(test_file_dir, 'reads-ds_2P.fastq'))]
    DB_KEY = 'virulencefinder'

    CONFIG_TEMPLATE = {
        'analyses': ['resfinder'],
        'detection_method': 'srst2',
        'gene_detection':
            {'virulencefinder':
                 {'filtering_method': 'cluster',
                  'min_coverage': 25,
                  'min_percent_identity': 40,
                  'path': '/data/srst2/gene_db/VirulenceFinder-Ecoli-clustered_80/',
                  'extra_column': ['Protein function', 'protein_function']}},
        'logging': False,
        'pipeline_job_id': 368,
        'pipeline_name': 'Test Pipeline',
        'sample_name': 'my_sample'
    }

    def setUp(self):
        """
        Sets up the resources before running the test.
        :return: None
        """
        self.running_dir = tempfile.mkdtemp(None, 'camel_', TestGeneDetection.camel.config['temp_dir'])

    def test_blastn(self):
        """
        Tests the gene detection using blastn.
        :return: None
        """
        # Get the target output file (HTML report)
        output_report_section = os.path.join(self.running_dir, OUTPUT_GENE_DETECTION_REPORT.format(
            db=TestGeneDetection.DB_KEY))

        # Copy input file
        input_file_new_path = os.path.join(self.running_dir, INPUT_GENE_DETECTION_FASTA.format(db="virulencefinder"))
        if not os.path.isdir(os.path.dirname(input_file_new_path)):
            os.makedirs(os.path.dirname(input_file_new_path))
        SnakemakeUtils.dump_object([TestGeneDetection.FILE_FASTA], input_file_new_path)

        # Create the config file
        config_path = os.path.join(self.running_dir, 'config.yml')
        with open(config_path, 'w') as handle:
            content = CONFIG_TEMPLATE.copy()
            content['detection_method'] = 'blast'
            content['working_dir'] = self.running_dir
            yaml.dump(content, handle)

        # Execute Snakemake
        SnakemakeTestUtils.run_snakemake(
            SNAKEFILE_GENE_DETECTION, config_path, self.running_dir, output_report_section, 8)

        # Save report
        report_path = os.path.join(self.running_dir, 'report.html')
        SnakemakeTestUtils.save_report(report_path, output_report_section)
        logging.info("Full path to report: {}".format(report_path))
        self.assertGreater(os.path.getsize(output_report_section), 0, "Output file is empty")

    def test_srst2(self):
        """
        Tests the gene detection using SRST2.
        :return: None
        """
        # Get the target output file (HTML report)
        output_report_section = os.path.join(self.running_dir, OUTPUT_GENE_DETECTION_REPORT.format(
            db=TestGeneDetection.DB_KEY))

        # Copy input file
        input_file_new_path = os.path.join(self.running_dir, INPUT_GENE_DETECTION_FASTQ_PE.format(db="virulencefinder"))
        if not os.path.isdir(os.path.dirname(input_file_new_path)):
            os.makedirs(os.path.dirname(input_file_new_path))
        SnakemakeUtils.dump_object(TestGeneDetection.FILES_FASTQ_PE, input_file_new_path)

        # Create the config file
        config_path = os.path.join(self.running_dir, 'config.yml')
        with open(config_path, 'w') as handle:
            content = CONFIG_TEMPLATE.copy()
            content['detection_method'] = 'srst2'
            content['working_dir'] = self.running_dir
            yaml.dump(content, handle)

        # Execute Snakemake
        SnakemakeTestUtils.run_snakemake(
            SNAKEFILE_GENE_DETECTION, config_path, self.running_dir, output_report_section, 8)

        # Save report
        report_path = os.path.join(self.running_dir, 'report.html')
        SnakemakeTestUtils.save_report(report_path, output_report_section)
        logging.info("Full path to report: {}".format(report_path))
        self.assertGreater(os.path.getsize(output_report_section), 0, "Output file is empty")

    def test_empty_report(self):
        """
        Tests the empty report generation.
        :return: None
        """
        # Get the target output file (HTML report)
        output_empty_report = os.path.join(self.running_dir, OUTPUT_GENE_DETECTION_REPORT_EMPTY.format(
            db=TestGeneDetection.DB_KEY))

        # Create the config file
        config_path = os.path.join(self.running_dir, 'config.yml')
        with open(config_path, 'w') as handle:
            content = CONFIG_TEMPLATE.copy()
            content['detection_method'] = 'blast'
            content['working_dir'] = self.running_dir
            yaml.dump(content, handle)

        # Execute Snakemake
        SnakemakeTestUtils.run_snakemake(
            SNAKEFILE_GENE_DETECTION, config_path, self.running_dir, output_empty_report, 8)

        # Save report
        report_path = os.path.join(self.running_dir, 'report.html')
        SnakemakeTestUtils.save_report(report_path, output_empty_report)
        logging.info("Full path to report: {}".format(report_path))
        self.assertGreater(os.path.getsize(output_empty_report), 0, "Output file is empty")


if __name__ == '__main__':
    unittest.main()
