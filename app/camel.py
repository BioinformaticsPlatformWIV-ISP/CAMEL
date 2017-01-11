from app.connection.connection import Connection
from app.loggers.logmanager import LogManager


class Camel(object):
    """
    Main class for camel.
    """

    def __init__(self, database_config):
        """
        Initializes a CAMEL system.
        """
        self._connection = Connection(database_config)
        LogManager.initialize()

    @property
    def connection(self):
        """
        Returns the database connection.
        :return: Database connection
        """
        return self._connection

