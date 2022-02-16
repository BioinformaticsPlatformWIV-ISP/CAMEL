import argparse
import logging
from pathlib import Path
from typing import Optional, Any, Dict, List, Tuple, Sequence, Union

import abc
import shutil

from camel.app.camel import Camel
from camel.app.components import mainscriptutils
from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.error.snakemakeexecutionerror import SnakemakeExecutionError
from camel.app.loggers import fileloggerutils
from camel.app.pipeline.pipeline import Pipeline
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils


class BasePipeline(object, metaclass=abc.ABCMeta):
    """
    This class is the base class for pipelines.
    """

    def __init__(self, name: str, version: str, snakefile: str, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the pipeline.
        :param name: Pipeline name
        :param version: Pipeline version
        """
        self._name = name
        self._version = version
        self._snakefile = snakefile
        self._args = self._parse_arguments(args)
        self._keep_logs = True if self._args.log else Camel.get_instance().config.get('logging', {}).get(
            'keep_logs', False)
        self._pipeline = Pipeline(name, Camel.get_instance(), self._args.log, self._args.log)

    @staticmethod
    @abc.abstractmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments. Should be implemented by the subclasses.
        :param args: Arguments (optional)
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
        argument_parser.add_argument('--fastq-pe', nargs=2, type=Path, help="FASTQ input files")
        argument_parser.add_argument('--fastq-pe-names', nargs=2, help="FASTQ input file names")

        # Output
        argument_parser.add_argument('--working-dir', type=Path, default=Path.cwd())

        # Options
        argument_parser.add_argument('--threads', default=8, type=int)
        argument_parser.add_argument(
            '--library', help="Adapter library that was used for the sequencing",
            choices=['NexteraPE', 'TruSeq2', 'TruSeq3'], default='NexteraPE')

        # Logging
        argument_parser.add_argument(
            '--galaxy-job-id', type=str, help='Job id of the run in galaxy (used for logging')
        argument_parser.add_argument(
            '--log', action='store_true', help="If this flag is set, config file and error logs are kept")

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

    @property
    def galaxy_job_id(self) -> Union[int, None]:
        """
        Returns the galaxy job id (if there is one).
        :return: Galaxy job id
        """
        return self._args.galaxy_job_id if 'galaxy_job_id' in self._args else None

    def _get_fastq_input_links(self) -> List[List[Tuple[Path, str]]]:
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
        dir_links = self._args.working_dir / 'input'
        if not dir_links.exists():
            dir_links.mkdir(parents=True)

        # Link files
        paths_new = []
        for path_orig, link_name in links:
            path_new = dir_links / link_name
            logging.debug(f"Symlinking input file: {path_orig} -> {link_name}")
            if path_new.is_symlink():
                path_new.unlink()
            path_new.symlink_to(path_orig)
            paths_new.append(path_new)

        # Return output dictionary
        return [{'name': p.name, 'path': p} for p in paths_new]

    def _run_snakemake_main(self, config_file: str) -> None:
        """
        Runs the main snakefile for the pipeline.
        :param config_file: Configuration file
        :return: None
        """
        # Store config file
        if self._keep_logs is True:
            fileloggerutils.store_config_file(Path(config_file), self._name, self.galaxy_job_id)

        # Clear existing Galaxy output files when html output is selected
        if 'output_html' in self._args:
            mainscriptutils.prepare_galaxy_output(Path(self._args.output_dir), Path(self._args.output_html))

        # Path to the logfile
        log_file = self._args.working_dir / 'camel.log'
        try:
            # Run snakemake
            SnakePipelineUtils.run_snakemake(
                self._snakefile, config_file, [], self._args.working_dir, self._args.threads)
            logging.info("Pipeline finished successfully")
        except SnakemakeExecutionError as err:
            if self._keep_logs and log_file.exists():
                log_file = fileloggerutils.store_log_file(log_file, self._name, self.galaxy_job_id, True)
                raise RuntimeError(f"Error executing Snakemake. Check log for more information: {log_file}")
            else:
                raise err

        # Copy log file to output directory if that directory is given
        if log_file.exists() and 'output_dir' in self._args:
            shutil.copyfile(str(log_file), str(Path(self._args.output_dir) / 'camel.log'))

        # Store log file
        if self._keep_logs and log_file.exists():
            fileloggerutils.store_log_file(log_file, self._name, self.galaxy_job_id)

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
            'working_dir': str(self._args.working_dir),
            'read_trimming': {'adapter': self._args.library}
        }
        return template_data
