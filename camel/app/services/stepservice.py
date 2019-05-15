from dataclasses import dataclass
from typing import Optional, Union

import snakemake
from psycopg2.extras import Json

from camel.app.connection.connection import Connection
from camel.app.services.service import Service


@dataclass(frozen=True)
class StepOutput:
    """
    This class is used to keep an output of a step.
    """
    pipeline_job_id: int
    rule_name: str
    type: str
    key: str
    index: int
    hash: str
    wildcards: Optional[snakemake.io.Wildcards] = None
    is_pipeline_io: bool = False

    @property
    def wildcards_json(self) -> Union[None, Json]:
        """
        Returns the wildcards as a JSON string.
        :return: Wildcards as string
        """
        if self.wildcards is None:
            return None
        return Json({k: v for k, v in self.wildcards.items()})


class StepService(Service):
    """
    This class will perform operations on the DB regarding steps.
    """

    def __init__(self, connection: Connection) -> None:
        """
        Initializes step service.
        :param connection: Connection to the database
        :return: None
        """
        super(StepService, self).__init__(connection)

    def log_output(self, step_output: StepOutput) -> None:
        """
        Logs the given output data for this step. The output data should be a list that contains (in order):
        pipeline job id, rule name, wildcards, type (e.g. file), key (e.g. FASTQ_PE), index of the file (i.e. first,
        second, ... file), hash of the file
        :param step_output: Step output
        :return: None
        """
        sql = """
        INSERT INTO logging.input_output(pipeline_job_id, rule_name, wildcards, type, key, index, hash, pipeline_io)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s);
        """
        self.db_connection.insert(sql, (
            step_output.pipeline_job_id, step_output.rule_name, step_output.wildcards_json, step_output.type,
            step_output.key, step_output.index, step_output.hash, step_output.is_pipeline_io))
