import datetime

from app.connection.connection import Connection
from app.services.service import Service


class PipelineService(Service):
    """
    This class will perform operations on the DB regarding pipelines. This service is attached to a single pipeline
    instance.
    """

    def __init__(self, pipeline_name: str, connection: Connection) -> None:
        """
        Initializes pipeline service.
        :param pipeline_name: Name of the pipeline
        :param connection: Connection to the database
        :return: None
        """
        super(PipelineService, self).__init__(connection)
        self._pipeline_id = self.__get_pipeline_id(pipeline_name)

    def __get_pipeline_id(self, name: str) -> int:
        """
        Returns the id of the pipeline with the given name.
        :param name: Pipeline name
        :return: Pipeline id
        """
        sql = """
        SELECT pipe_id
        FROM pipelines.pipeline
        WHERE name = %s;"""
        try:
            return self.db_connection.query(sql, (name,))[1][0]
        except IndexError:
            raise ValueError(f"Pipeline '{name}' not found in the database.")

    def insert_pipeline_job(self) -> int:
        """
        Inserts a new job into the database and returns the resulting id
        :return: Job id
        """
        date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")
        sql = """
          INSERT INTO pipelines.pipeline_job (date, pipeline_id)
          VALUES (%s, %s) RETURNING pipeline_job_id;
        """
        return self.db_connection.insert(sql, (date, self._pipeline_id))

    def log_initial_input(self, pipeline_job_id: int, type_: str, key: str, index: int, hash_value: str) -> None:
        """
        Logs the initial pipeline input.
        :return: None
        """
        sql = """
        INSERT INTO logging.input_output(pipeline_job_id, rule_name, wildcards, type, key, index, hash, pipeline_io)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s);
        """
        self.db_connection.insert(sql, (pipeline_job_id, 'initial', None, type_, key, index, hash_value, True))
