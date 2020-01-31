import argparse
import logging
from pathlib import Path
from typing import Optional, Any, Dict, List, Tuple

import abc
import os
import shutil

from camel.app.camel import Camel
from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.mainscripthelper import MainScriptHelper
from camel.app.error.snakemakeexecutionerror import SnakemakeExecutionError
from camel.app.pipeline.pipeline import Pipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils


class BasePipeline(object, metaclass=abc.ABCMeta):
    """
    This class is the base class for pipelines.
    """

    def __init__(self, name: str, version: str, snakefile: str, args: Optional[argparse.Namespace] = None) -> None:
        """
        Initializes the pipeline.
        :param name: Pipeline name
        :param version: Pipeline version
        """
        self._name = name
        self._version = version
        self._snakefile = snakefile
        self._args = self._parse_arguments() if args is None else args
        self._working_dir = Path(self._args.working_dir)
        self._pipeline = Pipeline(name, Camel.get_instance(), 'pipeline' if self._args.db_logging else None)

    @staticmethod
    @abc.abstractmethod
    def _parse_arguments() -> argparse.Namespace:
        """
        Parses the command line arguments. Should be implemented by the subclasses.
        :return: None
        """
        pass

    @staticmethod
    def add_common_arguments(argument_parser: argparse.ArgumentParser) -> None:
        """
        Adds the common arguments for the pipeline
        :param argument_parser: Argument parse
        :return: None
        """
        # Input
        argument_parser.add_argument('--sample-name', type=str)
        argument_parser.add_argument('--fastq-pe', nargs=2, help="FASTQ input files")
        argument_parser.add_argument('--fastq-pe-names', nargs=2, help="FASTQ input file names")

        # Output
        argument_parser.add_argument('--output-dir', required=True, type=str)
        argument_parser.add_argument('--output-html', required=True, type=str)
        argument_parser.add_argument('--output-tsv', help="Output file for the summary", required=True)
        argument_parser.add_argument('--working-dir', default=os.path.abspath('.'), type=str)

        # Options
        argument_parser.add_argument(
            '--detection-method', help="Type of allele detection: local alignment (blast), read mapping (srst2)",
            choices=['blast', 'srst2'], default='blast')
        argument_parser.add_argument('--threads', default=8, type=int)
        argument_parser.add_argument(
            '--report-include-fastq', help="Include the FASTQ files in the report", action='store_true')
        argument_parser.add_argument(
            '--report-include-bam', help="Include the BAM file in the report", action='store_true')
        argument_parser.add_argument(
            '--library', help="Adapter library that was used for the sequencing",
            choices=['nextera', 'truseq2', 'truseq3'], default='nextera')
        argument_parser.add_argument(
            '--db-logging', action='store_true', help="If this flag is set, output is logged to database")

    @property
    def name(self) -> str:
        """
        Returns the pipeline name.
        :return: Name
        """
        return self._name

    @property
    def title(self) -> str:
        """
        Returns the title of this pipeline, defaults to pipeline name.
        :return: Title
        """
        return self.name

    @property
    def version(self) -> str:
        """
        Returns the pipeline version.
        :return: Version
        """
        return self._version

    @property
    def sample_name(self) -> str:
        """
        Returns the sample name.
        :return: Sample name
        """
        if self._args.sample_name is not None:
            return FileSystemHelper.make_valid(self._args.sample_name)
        elif self._args.fastq_pe_names is not None:
            return FastqUtils.get_sample_name(self._args.fastq_pe_names[0])
        else:
            return FastqUtils.get_sample_name(self._args.fastq_pe[0])

    def _get_fastq_input_links(self) -> List[List[Tuple[str, str]]]:
        """
        Returns the links to the input FASTQ files.
        :return: Links
        """
        links = []
        for read_nb, path in enumerate(self._args.fastq_pe, start=1):
            gzipped = FileSystemHelper.is_gzipped(path)
            links.append([path, f"{self.sample_name}_{read_nb}.fastq{'.gz' if gzipped else ''}"])
        return links

    def _symlink_input(self) -> List[Dict[str, Any]]:
        """
        Symlinks the input files.
        :return: List of FASTQ input dictionaries
        """
        # Determine link names
        links = self._get_fastq_input_links()

        # Create directory
        dir_links = self._working_dir / 'input'
        if not dir_links.exists():
            dir_links.mkdir(parents=True)

        # Link files
        paths_new = []
        for path_orig, link_name in links:
            path_new = os.path.join(dir_links, link_name)
            logging.debug(f"Symlinking input file: {path_orig} -> {link_name}")
            if os.path.islink(path_new):
                os.remove(path_new)
            os.symlink(path_orig, path_new)
            paths_new.append(path_new)

        # Return output dictionary
        return [{'name': os.path.basename(p), 'path': p} for p in paths_new]

    def _run_snakemake_main(self, config_file: str) -> None:
        """
        Runs the main snakefile for the pipeline.
        :param config_file: Configuration file
        :return: None
        """
        self._pipeline.log_config_file(config_file)
        MainScriptHelper.prepare_galaxy_output(self._args.output_dir, self._args.output_html)
        try:
            SnakePipelineUtils.run_snakemake(
                self._snakefile, config_file, [], self._working_dir, self._args.threads)
            log_file = self._working_dir / 'camel.log'
            if log_file.exists():
                shutil.copyfile(str(log_file), str(Path(self._args.output_dir) / 'camel.log'))
            logging.info("Pipeline finished successfully")
        except SnakemakeExecutionError as err:
            self._pipeline.log_error_to_file(err)

    def get_template_data(self, input_key: str, input_data: [List[Dict[str, str]]]) -> Dict[str, Any]:
        """
        Returns the template data that is common to all pipeline.
        :param input_key: FASTQ input key
        :param input_data: FASTQ input files
        :return: Template data
        """
        template_data = {
            'pipeline': {
                'name': self._name,
                'version': f"{self._version}",
                'title': self.title
            },
            input_key: input_data,
            'sample_name': self.sample_name,
            'output_report': self._args.output_html,
            'output_tabular': self._args.output_tsv,
            'output_dir': self._args.output_dir,
            'working_dir': str(self._working_dir),
            'detection_method': self._args.detection_method,
            'read_trimming': {'export_fastq': self._args.report_include_fastq}
        }
        if self._pipeline.is_logged:
            template_data['logging_level'] = self._pipeline.logging_level
            template_data['pipeline_job_id'] = self._pipeline.job_id
        return template_data
