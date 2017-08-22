import logging

from app.pipeline.pipeline import Pipeline
from app.services.pipelineservice import PipelineService


class SnakePipeline(Pipeline):
    """
    Class meant to handle the workflow of steps.
    """

    def __init__(self, name, camel, db_logging=False, job_id=None):
        """
        Initializes a pipeline.
        :param name: Name of the pipeline
        :param camel: CAMEL instance
        :param db_logging: If True, inputs & outputs are logged in the database.
        """
        super(SnakePipeline, self).__init__(camel, db_logging)
        self._name = name
        if db_logging is True:
            self._pipeline_service = PipelineService(self._name, camel.connection)
            if job_id is None:
                self._job_id = self._pipeline_service.insert_pipeline_job()
            else:
                self._job_id = job_id
        logging.info("Created pipeline {}".format(self._name))

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

    def set_initial_input(self, files):
        """
        Sets the initial inputs for the pipeline allowing it to be logged if logging is requested.
        :param files: Dictionary of input files
        :return: None
        """
        super(SnakePipeline, self).set_initial_input(files)
        if self.db_logging is True:
            self._log_initial_input()

    def get_initial_input(self, key=None):
        """
        Returns a list of file name strings so that they can be used in Snakemake. If a key is given only the files
        for that key are returned, otherwise all files in the dictionary are returned.
        :param key: Optional input file key
        :return: List of input files
        """
        iofiles = []
        if key is None:
            for values in self._initial_input.values():
                iofiles += values
        else:
            iofiles = self._initial_input[key]
        files = []
        for file_ in iofiles:
            files.append(file_.path)
        return files
