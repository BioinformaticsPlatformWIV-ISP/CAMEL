import logging

from abc import ABC, abstractmethod


class Pipeline(ABC):
    """
    Class meant to handle the workflow of steps.
    """

    def __init__(self, camel, db_logging=False):
        """
        Initializes a pipeline.
        :param camel: CAMEL instance
        :param db_logging: If True, inputs & outputs are logged in the database.
        """
        self._camel = camel
        self._db_logging = db_logging
        self._name = None
        self._initial_input = None
        self._configs = None
        self._pipeline_service = None
        self._job_id = None

    @property
    def job_id(self):
        """
        Returns the job id of this pipeline.
        :return: Job id
        """
        return self._job_id

    @property
    def pipeline_service(self):
        """
        Returns the pipeline service
        :return: Pipeline service
        """
        return self._pipeline_service

    def set_configs(self, configs):
        """
        Sets up the configuration of the pipeline.
        :param configs: Configuration
        :return: None
        """
        logging.info("Pipeline configuration: {}".format(configs))
        self._configs = configs

    def _log_initial_input(self):
        """
        Logs the initial input of the pipeline.
        :return: None
        """
        for key, files in self._initial_input.items():
            for i in range(0, len(files)):
                if files[i].logged:
                    self._pipeline_service.log_initial_input(self._job_id, files[i].TYPE_NAME, key, i, files[i].hash)
                    logging.debug('Initial input {} ({}) logged'.format(key, i))

    @abstractmethod
    def set_initial_input(self, files):
        """
        Sets the initial input files of the pipeline.
        :param files: dictionary of files to import
        :return: None
        """
        if type(files) is dict:
            self._initial_input = files
        else:
            raise TypeError("Input object should be a dictionary")
