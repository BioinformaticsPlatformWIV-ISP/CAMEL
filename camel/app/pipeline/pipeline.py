import gzip
import logging
from itertools import chain
from pathlib import Path
from typing import Optional

import shutil

from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.error.snakemakeexecutionerror import SnakemakeExecutionError


class Pipeline(object):
    """
    Class meant to handle the workflow of steps.
    """

    def __init__(self, name: str, camel: Camel, keep_config: bool = False, keep_error_log: bool = False) -> None:
        """
        Initializes a pipeline.
        :param name: Pipeline name
        :param camel: CAMEL instance
        :param keep_config: If True, pipeline config is saved on disk
        :param keep_error_log: If True, log is saved in case of error
        """
        self._camel = camel
        self._name = name
        self._initial_input = None
        self._job_id = None
        self._keep_config = keep_config
        self._keep_error_log = keep_error_log
        logging.info(f"Initialized pipeline '{self._name}'")

    @property
    def name(self) -> str:
        """
        Returns the pipeline name.
        :return: Name
        """
        return self._name

    @property
    def keep_config(self) -> bool:
        """
        Returns True if the config file needs to be kept, False otherwise.
        :return: True if config is kept
        """
        return self._keep_config

    @property
    def keep_error_log(self) -> bool:
        """
        Returns True if the error log needs to be kept.
        :return: True if error log is kept
        """
        return self._keep_error_log

    def set_initial_input(self, files: dict) -> None:
        """
        Sets the initial inputs for the pipeline allowing it to be logged if logging is requested.
        :param files: Dictionary of input files
        :return: None
        """
        self._initial_input = files
        self._log_initial_input()

    def get_initial_input(self, key: str = None) -> list:
        """
        Returns a list of file paths so that they can be used in Snakemake. If a key is given only the file paths
        for that key are returned, otherwise all file paths in the dictionary are returned.
        :param key: Optional input file key
        :return: List of input file paths
        """
        io_files = list(chain(*self._initial_input.values())) if key is None else self._initial_input[key]
        return [file_.path for file_ in io_files]

    def _log_initial_input(self) -> None:
        """
        Logs the initial input of the pipeline.
        :return: None
        """
        for key, files in self._initial_input.items():
            for index, io in enumerate(files):
                if not io.is_logged:
                    continue
                logging.info(f"Pipeline input log: type={io.type_name} key={key} index={index} hash={io.hash}")

    def log_config_file(self, config_file: str, galaxy_job_id: Optional[str] = None) -> None:
        """
        Exports the config file from the pipeline.
        :param galaxy_job_id: Galaxy job id
        :param config_file: Config file to export
        :return: None
        """
        export_path = Path(Camel.get_instance().config['config_dump_dir']) / '{}.yml.gz'.format('_'.join([
            FileSystemHelper.make_valid(self.name),
            FileSystemHelper.get_timestamp_str(),
            galaxy_job_id if galaxy_job_id is not None else 'NA',
        ]))
        with gzip.open(export_path, 'wb') as file_out, open(config_file, 'rb') as file_in:
            shutil.copyfileobj(file_in, file_out)
        logging.info(f"Config file exported to '{export_path}'")

    def log_error_to_file(self, error: SnakemakeExecutionError) -> None:
        """
        Dumps the pipeline error log.
        :param error: Error raised by Snakemake
        :return: None
        """
        output_logfile = Path(Camel.get_instance().config['error_log_dir']) / 'error-{}-{}.txt'.format(
            FileSystemHelper.make_valid(self.name).lower(), FileSystemHelper.get_timestamp_str())
        logging.debug(f"Dumping error log: {output_logfile}")
        with open(output_logfile, 'w') as handle_in:
            handle_in.write('Stdout:\n')
            handle_in.write(error.stdout)
            handle_in.write('-' * 10 + '\n')
            handle_in.write('Stderr:\n')
            handle_in.write(error.stderr)
        raise RuntimeError(f"Error executing Snakemake. Check log for more information: {output_logfile}")
