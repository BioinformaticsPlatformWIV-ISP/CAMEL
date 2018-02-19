import logging

import os
import yaml

from app.components.filesystemhelper import FileSystemHelper
from app.connection.connection import Connection
from app.loggers.logmanager import LogManager
from app.services.basetoolservice import BaseToolService
from app.services.dbtoolservice import DbToolService
from app.services.yamltoolservice import YAMLToolService
from config import DB_CONFIG, LOGGING_CONFIG, MAIN_CONFIG
from tool_data import TOOL_DATA_DIR


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

    @staticmethod
    def get_tool_data_path(tool_name: str, tool_version: str) -> str:
        """
        Returns the path of the tool data for the tool with the given name and version.
        :param tool_name: Tool name
        :param tool_version: Tool version
        :return: Path
        """
        return os.path.join(TOOL_DATA_DIR, '{}-{}.yml'.format(
            FileSystemHelper.make_valid(tool_name).lower(),
            FileSystemHelper.make_valid(tool_version)))

    def get_tool_service(self, tool_name: str, tool_version: str) -> BaseToolService:
        """
        Returns the tool service for the tool with the given name and version.
        :return: Tool service
        """
        logging.debug('Retrieving tool service')
        source = self._config.get('tool_service', 'db')
        if source == 'db':
            return DbToolService(tool_name, tool_version, self.connection)
        elif source == 'yaml':
            tool_data_path = Camel.get_tool_data_path(tool_name, tool_version)
            if not os.path.isfile(tool_data_path):
                raise FileNotFoundError('Tool data file not found: {}'.format(os.path.basename(tool_data_path)))
            return YAMLToolService(tool_data_path)
        else:
            raise ValueError("Invalid 'tool_service' value: {}".format(source))
