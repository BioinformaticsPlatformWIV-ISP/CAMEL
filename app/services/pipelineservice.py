import logging

import datetime

from app.parameter.parameter import Parameter
from app.services.service import Service


class PipelineService(Service):
    """
    This class will perform operations on the DB regarding pipelines. This service is attached to a single pipeline
    instance.
    """

    def __init__(self, pipeline_name, connection):
        """
        Initializes pipeline service.
        :param pipeline_name: Name of the pipeline
        :param connection: Connection to the database
        :return: None
        """
        super(PipelineService, self).__init__(connection)
        self._pipeline_id = self.__get_pipeline_id(pipeline_name)

    def __get_pipeline_id(self, name):
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
            return self.db_connection.query(sql, [name])[1][0]
        except IndexError:
            raise ValueError("Pipeline '{}' not found in the database.".format(name))

    def insert_pipeline_job(self):
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

    def insert_step(self, name):
        """
        Inserts a step in to the database.
        :param name: Step name
        :return: Step id
        """
        sql = """
        INSERT INTO pipelines.pipeline_step(name, pipe_id)
        VALUES(%s, %s)
        RETURNING pipeline_step_id;"""
        logging.info("Inserting pipeline step '{}' in the database".format(name))
        return self.db_connection.insert(sql, [name, self._pipeline_id])

    def get_step_ids(self):
        """
        Returns the ids of the pipeline steps in the database.
        :return: Dictionary containing step name and step id.
        """
        sql = """
        SELECT name, pipeline_step_id
        FROM pipelines.pipeline_step
        WHERE pipe_id = %s;"""
        return {result[0]: result[1] for result in self.db_connection.query(sql, [self._pipeline_id])[1:]}

    def get_step_id(self, name):
        """
        Returns the id of the pipeline step with the given name in the database.
        :param name: Name of the step
        :return: Step id if the step exists, else None
        """
        sql = """
        SELECT pipeline_step_id
        FROM pipelines.pipeline_step
        WHERE pipe_id = %s
        AND name = %s;"""
        results = self.db_connection.query(sql, [self._pipeline_id, name])
        return results[1][0] if len(results) > 1 else None

    def get_pipeline_parameters(self):
        """
        Returns the pipeline parameters.
        :return: Pipeline parameters list
        """
        sql = """
        SELECT s.name, p.name, p.option, sp.value, sp.disabled
        FROM pipelines.step_tools_parameter sp, tools.tool_parameter p, pipelines.pipeline_step s
        WHERE sp.tool_parameter_id = p.tool_parameter_id
        AND sp.pipeline_step_id IN (
          SELECT pipeline_step_id FROM pipelines.pipeline_step
          WHERE pipe_id = %s)
        AND sp.active is True
        AND sp.pipeline_step_id = s.pipeline_step_id
        ORDER BY sp.p_index;
        """
        parameters = []
        for step_name, name, option, value, disabled in self.db_connection.query(sql, [self._pipeline_id])[1:]:
            if disabled:
                value = False
            parameters.append([step_name, Parameter(name, option, value)])
        return parameters

    def get_pipeline_step_parameters(self, name):
        """
        Returns the pipeline parameters for step with the specified name.
        :param name: Name of the step
        :return: Pipeline parameters list
        """
        sql = """
        SELECT s.name, p.name, p.option, sp.value, sp.disabled
        FROM pipelines.step_tools_parameter sp, tools.tool_parameter p, pipelines.pipeline_step s
        WHERE sp.tool_parameter_id = p.tool_parameter_id
        AND sp.pipeline_step_id IN (
          SELECT pipeline_step_id FROM pipelines.pipeline_step
          WHERE pipe_id = %s
          AND name = %s)
        AND sp.active is True
        AND sp.pipeline_step_id = s.pipeline_step_id
        ORDER BY sp.p_index;
        """
        parameters = []
        for step_name, name, option, value, disabled in self.db_connection.query(sql, [self._pipeline_id, name])[1:]:
            if disabled:
                value = False
            parameters.append(Parameter(name, option, value))
        return parameters

    def log_initial_input(self, pipeline_job_id, type_, key, index, hash_value):
        """
        Logs the initial pipeline input.
        :return: None
        """
        sql = """
        INSERT INTO logging.pipeline_input(pipeline_job_id, type, key, index, hash)
        VALUES (%s, %s, %s, %s, %s);
        """
        self.db_connection.insert(sql, [pipeline_job_id, type_, key, index, hash_value])
