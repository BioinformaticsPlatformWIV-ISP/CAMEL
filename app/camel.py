from app.connection.connection import Connection
from app.loggers.logmanager import LogManager
from config import DB_CONFIG, LOGGING_CONFIG


class Camel(object):
    """
    Main class for camel.
    """

    def __init__(self, database_config=DB_CONFIG, logging_config=LOGGING_CONFIG):
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
