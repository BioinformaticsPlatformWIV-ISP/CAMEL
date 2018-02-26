import logging

from camel.app.camel import Camel
from camel.app.services.pipelineservice import PipelineService


class Pipeline(object):
    """
    Class meant to handle the workflow of steps.
    """

    def __init__(self, name: str, camel: Camel, db_logging: bool=False, job_id: int=None) -> None:
        """
        Initializes a pipeline.
        :param camel: CAMEL instance
        :param db_logging: If True, inputs & outputs are logged in the database.
        """
        self._camel = camel
        self._name = name
        self._initial_input = None
        self._configs = None
        self._job_id = None
        self._db_logging = db_logging
        if db_logging:
            self._pipeline_service = PipelineService(self._name, camel.connection)
            self._job_id = self._pipeline_service.insert_pipeline_job() if job_id is None else job_id
        else:
            self._pipeline_service = None
        logging.info("Created pipeline {}".format(self._name))

    @property
    def job_id(self):
        """
        Returns the job id of this pipeline.
        :return: Job id
        """
        return self._job_id

    @property
    def name(self):
        """
        Returns the pipeline name.
        :return: Name
        """
        return self._name

    @property
    def pipeline_service(self):
        """
        Returns the pipeline service
        :return: Pipeline service
        """
        return self._pipeline_service

    @property
    def configs(self):
        """
        Returns the pipeline configs
        :return: Pipeline configs
        """
        return self._configs

    @property
    def db_logging(self):
        """
        Returns the boolean that says whether or not to log in the database
        :return: True/False
        """
        return self._db_logging

    def set_initial_input(self, files: dict) -> None:
        """
        Sets the initial inputs for the pipeline allowing it to be logged if logging is requested.
        :param files: Dictionary of input files
        :return: None
        """
        self._initial_input = files
        if self.db_logging:
            self._log_initial_input()

    def get_initial_input(self, key: str=None) -> list:
        """
        Returns a list of file name strings so that they can be used in Snakemake. If a key is given only the files
        for that key are returned, otherwise all files in the dictionary are returned.
        :param key: Optional input file key
        :return: List of input files
        """
        io_files = []
        if key is None:
            for values in self._initial_input.values():
                io_files += values
        else:
            io_files = self._initial_input[key]
        files = []
        for file_ in io_files:
            files.append(file_.path)
        return files

    def set_configs(self, configs: dict) -> None:
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
                    self._pipeline_service.log_initial_input(self._job_id, files[i].type_name, key, i, files[i].hash)
                    logging.debug('Initial input {} ({}) logged'.format(key, i))
