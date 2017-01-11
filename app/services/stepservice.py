from app.services.service import Service


class StepService(Service):
    """
    This class will perform operations on the DB regarding steps. This service is attached to a single step instance.
    """

    def __init__(self, step_id, connection):
        """
        Initializes step service.
        :param step_id: Step id
        :param connection: Connection to the database
        :return: None
        """
        super(StepService, self).__init__(connection)
        self._step_id = step_id

    def log_output(self, output_data):
        """
        Logs the given output data for this step.
        :param output_data: Output data
        :return: None
        """
        sql = """
        INSERT INTO logging.input_output(pipeline_job_id, pipeline_step_id, type, key, index, hash)
        VALUES (%s, %s, %s, %s, %s, %s);
        """
        self.db_connection.insert(sql, output_data)

    def log_job_parameter(self, parameter_id, pipeline_step_id, pipeline_job_id, value):
        """
        Logs a job parameter to the database.
        :return: None
        """
        sql = """
        INSERT INTO pipelines.job_step_tools_parameter (parameter_id, pipeline_step_id, pipeline_job_id, value,
          disabled)
        VALUES (%s, %s, %s, %s, %s)"""
        self.db_connection.insert(sql, [parameter_id, pipeline_step_id, pipeline_job_id, value, value is False])

    def get_parameter_id(self, tool_id, parameter_name):
        sql = """
        SELECT tp.tool_parameter_id
        FROM  tools.tool_parameter tp, tools.tool t
        WHERE tp.tool_id = t.tool_id
        AND t.tool_id = %s
        AND tp.name = %s;"""
        return self.db_connection.query(sql, [tool_id, parameter_name])[1][0]
