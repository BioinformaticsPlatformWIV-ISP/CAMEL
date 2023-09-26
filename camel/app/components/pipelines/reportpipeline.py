from pathlib import Path
from typing import Dict, Any, List

import abc
import argparse

from camel.app.components import mainscriptutils
from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.files.fileutils import FileUtils
from camel.app.components.phylogeny.snpphylogenyutils import InvalidInputError
from camel.app.components.pipelines.basepipeline import BasePipeline
from camel.app.loggers import logger


class ReportPipeline(BasePipeline, metaclass=abc.ABCMeta):
    """
    Baseclass for pipelines with a report (HTML) and tabular (TSV) output.
    """

    @staticmethod
    def add_common_arguments(argument_parser: argparse.ArgumentParser) -> None:
        """
        Adds the common arguments for the pipeline.
        :param argument_parser: Argument parser
        :return: None
        """
        BasePipeline.add_common_arguments(argument_parser)

        # Output
        argument_parser.add_argument('--output-dir', required=True, type=Path)
        argument_parser.add_argument('--output-html', required=True, type=Path)
        argument_parser.add_argument('--output-tsv', help="Output file for the summary", required=True, type=Path)

        # Options
        argument_parser.add_argument(
            '--detection-method', help="Type of allele detection: local alignment (blast), read mapping (srst2)",
            choices=['blast', 'kma', 'srst2'], default='blast')
        argument_parser.add_argument(
            '--report-include-fastq', help="Include the FASTQ files in the report", action='store_true')
        argument_parser.add_argument(
            '--report-include-bam', help="Include the BAM file in the report", action='store_true')

        # Parameters
        argument_parser.add_argument(
            '--cov-max', default=100.0, type=float,
            help='Maximum coverage (datasets with higher estimated coverage will be downsampled to the given value)')

    def get_template_data(self, input_key: str, input_data: [List[Dict[str, str]]]) -> Dict[str, Any]:
        """
        Returns the template data that is common to all pipeline.
        :param input_key: FASTQ input key
        :param input_data: FASTQ input files
        :return: Template data
        """
        # Convert Path objects to string in the config file
        config_data = super().get_template_data(input_key, [{
            k: str(v) if isinstance(v, Path) else v for k, v in fq_entry.items()} for fq_entry in input_data])
        mainscriptutils.dict_merge(config_data, {
            'output_dir': str(self._args.output_dir),
            'output_report': str(self._args.output_html),
            'output_tabular': str(self._args.output_tsv),
            'read_type': self._args.read_type,
            'detection_method': self._args.detection_method,
            'read_trimming': {'export_fastq': self._args.report_include_fastq}
        })
        if self._args.library is not None:
            config_data['read_trimming']['adapter'] = self._args.library
        return config_data

    def _validate_fastq_input(self) -> None:
        """
        Checks if the provided FASTQ input is valid.
        :return: None
        """
        logger.info(f'Checking FASTQ input')
        if self._args.read_type == 'illumina':
            nb_reads_fwd = FastqUtils.count_reads(self._args.fastq_pe[0])
            nb_reads_rev = FastqUtils.count_reads(self._args.fastq_pe[1])
            if not nb_reads_fwd == nb_reads_rev:
                raise InvalidInputError(
                    f'The number of forward ({nb_reads_fwd:,}) and reverse ({nb_reads_rev:,}) reads should be equal, '
                    'check that the input files provided are complete and correctly paired.')
            logger.info(f'FASTQ input is valid')
            logger.info(f'PE forward FASTQ hash: {FileUtils.hash_file(self._args.fastq_pe[0])}')
            logger.info(f'PE reverse FASTQ hash: {FileUtils.hash_file(self._args.fastq_pe[1])}')
        else:
            logger.debug('FASTQ checking not implemented yet')
