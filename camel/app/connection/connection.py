from typing import List, Optional, Tuple, Any

import psycopg2
import yaml


class Connection(object):
    """
    Class meant to create a connection to the database.
    """

    def __init__(self, config_path: str, db='local') -> None:
        """
        Initializes the connection object with user credentials
        Database configuration is read from the users folder, from a YAML file
        Cross OS compatible due to expanduser function
        :param config_path: Location of the configuration file with the connection parameters
        :param db: Database to connect to
        """
        with open(config_path) as f:
            self.config = yaml.safe_load(f)
        self._db = db

    @property
    def host(self) -> str:
        """
        Returns the host.
        :return: Host
        """
        return self.config[self._db]['host']

    @property
    def db_name(self) -> str:
        """
        Returns the database name.
        :return: Database name
        """
        return self.config[self._db]['dbname']

    def get_connection(self):
        """
        Returns a database connection.
        :return: Connection
        """
        conn = psycopg2.connect(host=self.config[self._db]['host'],
                                dbname=self.config[self._db]['dbname'],
                                user=self.config[self._db]['user'],
                                password=self.config[self._db]['password'])
        conn.autocommit = True
        return conn

    def query(self, query: str, params: Optional[Tuple[Any]] = None) -> List:
        """
        Method that initializes a db cursor and calls execute to perform the query
        :param query: sql query
        :param params: optional parameters to send in the query
        :return: Query result
        """
        with self.get_connection() as conn:
            cursor = conn.cursor()
            self.execute(cursor, query, params)
            return [[x[0].title() for x in cursor.description]] + [r for r in cursor.fetchall()]

    def insert(self, query: str, params: Optional[Tuple[Any]] = None) -> Optional[Any]:
        """
        Method that initializes a db cursor and calls execute to perform an insert query
        :param query: sql query
        :param params: parameters to send in the query
        :return: Returned db value or None depending on the query
        """
        with self.get_connection() as conn:
            cursor = conn.cursor()
            self.execute(cursor, query, params)
            return cursor.fetchone()[0] if 'RETURNING' in query else None

    @staticmethod
    def execute(cursor, query: str, params: Optional[Tuple[Any]] = None) -> None:
        """
        Executes sql code with proper exception handling
        :param cursor: db cursor just created
        :param query: sql code
        :param params: optional parameters to send in the query
        :return:
        """
        try:
            cursor.execute(query, params)
        except psycopg2.Error as e:
            raise ValueError("Error executing SQL query: {}".format(e))
