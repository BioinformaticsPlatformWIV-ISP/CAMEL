import os

import logging

from app.connection.connection import Connection


class Camel(object):
    """
    Main class for camel.
    """

    def __init__(self, database_config):
        """
        Initializes a CAMEL system.
        """
        self._connection = Connection(database_config)
        self._setup_logger()

    @property
    def connection(self):
        """
        Returns the database connection.
        :return: Database connection
        """
        return self._connection

    def _setup_logger(self):
        """
        Sets up the logger for this CAMEL instance.
        :return: None
        """
        formatter = logging.Formatter('CAMEL - %(levelname)s - %(asctime)s - %(message)s')
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        console_handler.setLevel(logging.DEBUG)
        logging.getLogger().addHandler(console_handler)
        file_handler = logging.FileHandler(os.path.join(os.getcwd(), 'camel.log'))
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logging.getLogger().addHandler(file_handler)
        logging.getLogger().setLevel(logging.DEBUG)

