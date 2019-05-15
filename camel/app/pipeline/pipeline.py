import gzip
import logging
from itertools import chain
from typing import Optional

import os
import shutil

from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.error.snakemakeexecutionerror import SnakemakeExecutionError
from camel.app.services.pipelineservice import PipelineService


class Pipeline(object):
    """
    Class meant to handle the workflow of steps.
    """

    def __init__(self, name: str, camel: Camel, logging_level: str = None, job_id: int = None) -> None:
        """
        Initializes a pipeline.
        :param camel: CAMEL instance
        :param logging_level: Logging level ('step', 'pipeline' or None)
        """
        self._camel = camel
        self._name = name
        self._initial_input = None
        self._job_id = None
        if not self.is_valid_logging_level(logging_level):
            raise ValueError(f"Logging level '{logging_level}' is not a valid logging level!")
        self._logging_level = logging_level
        self._db_logging = logging_level in ('pipeline', 'step')
        if self._db_logging:
            self._pipeline_service = PipelineService(self._name, camel.connection)
            self._job_id = self._pipeline_service.insert_pipeline_job() if job_id is None else job_id
        else:
            self._pipeline_service = None
        logging.info("Created pipeline '{}'".format(self._name))

    @property
    def job_id(self) -> Optional[int]:
        """
        Returns the job id of this pipeline.
        :return: Job id
        """
        return self._job_id

    @property
    def name(self) -> str:
        """
        Returns the pipeline name.
        :return: Name
        """
        return self._name

    @property
    def pipeline_service(self) -> Optional[PipelineService]:
        """
        Returns the pipeline service
        :return: Pipeline service
        """
        return self._pipeline_service

    @property
    def logging_level(self) -> str:
        """
        Returns the logging level of this pipeline.
        :return: Logging level
        """
        return self._logging_level

    @staticmethod
    def is_valid_logging_level(logging_level: str) -> bool:
        """
        Checks whether the given logging level is valid
        :param logging_level: Logging level to be checked
        :return: True if the level is allowed
        """
        return logging_level in ('pipeline', 'step', None)

    def set_initial_input(self, files: dict) -> None:
        """
        Sets the initial inputs for the pipeline allowing it to be logged if logging is requested.
        :param files: Dictionary of input files
        :return: None
        """
        self._initial_input = files
        if self._db_logging:
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
            for i in range(0, len(files)):
                if files[i].is_logged:
                    self._pipeline_service.log_initial_input(self._job_id, files[i].type_name, key, i, files[i].hash)
                    logging.debug('Initial input {} ({}) logged'.format(key, i))

    def log_config_file(self, config_file: str, galaxy_job_id: Optional[str] = None) -> None:
        """
        Exports the config file from the pipeline.
        :param galaxy_job_id: Galaxy job id
        :param config_file: Config file to export
        :return: None
        """
        if self._db_logging is False:
            logging.info("Logging disabled, config file not exported")
            return
        export_path = os.path.join(Camel.get_instance().config['config_dump_dir'], '{}.yml.gz'.format('_'.join([
            FileSystemHelper.get_timestamp_str(),
            galaxy_job_id if galaxy_job_id is not None else 'NA',
            str(self.job_id) if self.job_id is not None else 'NA',
            FileSystemHelper.make_valid(self.name)
        ])))
        with gzip.open(export_path, 'wb') as file_out, open(config_file, 'rb') as file_in:
            shutil.copyfileobj(file_in, file_out)
        logging.info(f"Config file exported to '{export_path}'")

    def log_error_to_file(self, error: SnakemakeExecutionError) -> None:
        """
        Dumps the pipeline error log.
        :param error: Error raised by Snakemake
        :return: None
        """
        if self._db_logging is False:
            logging.info("Logging disabled, error log not exported")
            output_logfile = 'NA (logging to file is disabled)'
        else:
            logging.debug("Dumping error log")
            output_logfile = os.path.join(Camel.get_instance().config['error_log_dir'], 'error-{}-{}.txt'.format(
                FileSystemHelper.make_valid(self.name).lower(), FileSystemHelper.get_timestamp_str()))
            with open(output_logfile, 'w') as handle_in:
                handle_in.write('Stdout:\n')
                handle_in.write(error.stdout)
                handle_in.write('-' * 10 + '\n')
                handle_in.write('Stderr:\n')
                handle_in.write(error.stderr)
        raise RuntimeError(f"Error executing Snakemake. Check log for more information: {output_logfile}")
