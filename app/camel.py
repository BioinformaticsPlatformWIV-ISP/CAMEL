from app.connection.connection import Connection
from app.loggers.logmanager import LogManager


class Camel(object):
    """
    Main class for camel.
    """

    def __init__(self, database_config, logging_config):
        """
        Initializes a CAMEL system.
        """
        LogManager.initialize(logging_config)
        self._connection = Connection(database_config)

    @property
    def connection(self):
        """
        Returns the database connection.
        :return: Database connection
        """
        return self._connection

