import yaml

from camel.app.connection.connection import Connection
from camel.app.loggers.logmanager import LogManager
from camel.config import DB_CONFIG, LOGGING_CONFIG, MAIN_CONFIG


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

        with open(MAIN_CONFIG) as f:
            self._config = yaml.safe_load(f)

    @property
    def connection(self):
        """
        Returns the database connection.
        :return: Database connection
        """
        return self._connection

    @property
    def config(self):
        """
        Returns the main config as specified in app/config/main.yml
        :return: Dict
        """
        return self._config
