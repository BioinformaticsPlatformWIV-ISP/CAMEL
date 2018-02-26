from collections import OrderedDict
from typing import Dict, Tuple, Optional

from psycopg2.extras import Json
from snakemake.io import Wildcards

from camel.app.connection.connection import Connection
from camel.app.services.service import Service


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

    def log_output(self, output_data: Tuple[int, str, Wildcards, str, str, int, str, bool]) -> None:
        """
        Logs the given output data for this step. The output data should be a list that contains (in order):
        pipeline job id, rule name, wildcards, type (e.g. file), key (e.g. FASTQ_PE), index of the file (i.e. first,
        second, ... file), hash of the file
        :param output_data: Output data
        :return: None
        """
        if output_data[2] is not None:
            output_data = output_data[:2] + (Json(self.__wildcards_to_dict(output_data[2])),) + output_data[3:]
        sql = """
        INSERT INTO logging.input_output(pipeline_job_id, rule_name, wildcards, type, key, index, hash, pipeline_io)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s);
        """
        self.db_connection.insert(sql, output_data)

    @staticmethod
    def __wildcards_to_dict(wildcards: Wildcards) -> Optional[Dict[str, str]]:
        """
        Extracts the wildcard names and values from the Snakemake Wildcards object and returns them as a dictionary.
        :param wildcards: Wildcards object from snakemake
        :return: Dictionary with wildcards names (key) and values (value)
        """
        if wildcards is None:
            return None
        dict_ = OrderedDict()
        for key, value in sorted(wildcards.items()):
            dict_[key] = value
        return dict_
