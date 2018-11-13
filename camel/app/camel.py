from typing import Optional

import yaml

from camel.app.connection.connection import Connection
from camel.app.loggers.logmanager import LogManager
from camel.config import DB_CONFIG, LOGGING_CONFIG, MAIN_CONFIG


class Camel(object):
    """
    Main class for camel.
    """

    _current_instance = None
    _logger_is_initialized = False

    def __init__(self, database_config: str=DB_CONFIG, logging_config: Optional[str]=LOGGING_CONFIG,
                 tool_parameter_loc: str=None) -> None:
        """
        Initializes a CAMEL system.
        :param database_config: Location of database config file
        :param logging_config: Location of logging config file
        :param tool_parameter_loc: Location of tool parameter YAML files
        """
        if not Camel._logger_is_initialized and logging_config is not None:
            LogManager.initialize(logging_config)
            Camel._logger_is_initialized = True
        self._connection = Connection(database_config)

        with open(MAIN_CONFIG) as f:
            self._config = yaml.safe_load(f)

        if self._config.get('tool_service', 'db') == 'yaml' and 'tool_parameter_loc' not in self._config:
            self._config['tool_parameter_loc'] = tool_parameter_loc

    @property
    def connection(self) -> Connection:
        """
        Returns the database connection.
        :return: Database connection
        """
        return self._connection

    @property
    def config(self) -> dict:
        """
        Returns the main config as specified in app/config/main.yml
        :return: Dict
        """
        return self._config

    @staticmethod
    def get_instance() -> 'Camel':
        """
        This method can be used to avoid recreating the CAMEL object multiple times.
        :return: Initialized CAMEL instance
        """
        if Camel._current_instance is None:
            Camel._current_instance = Camel()
        return Camel._current_instance
