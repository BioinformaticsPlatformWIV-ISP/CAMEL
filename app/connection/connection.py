import logging
import psycopg2
import yaml


class Connection(object):
    """
    Class meant to create a connection to the database.
    """

    def __init__(self, config, db='local'):
        """
        Initializes the connection object with user credentials
        Database configuration is read from the users folder, from a YAML file
        Cross OS compatible due to expanduser function
        :param config: Location of the configuration file with the connection parameters
        :param db: Database to connect to
        """
        with open(config) as f:
            self.config = yaml.safe_load(f)
        self._db = db

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

    def query(self, query, params=None):
        """
        Method that initializes a db cursor and calls execute to perform the query
        :param query: sql query
        :param params: optional parameters to send in the query
        :return:
        """
        with self.get_connection() as conn:
            cursor = conn.cursor()
            self.execute(cursor, query, params)
            return [[x[0].title() for x in cursor.description]] + [r for r in cursor.fetchall()]

    def insert(self, query, params):
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
    def execute(cursor, query, params=None):
        """
        Executes sql code with proper exception handling
        :param cursor: db cursor just created
        :param query: sql code
        :param params: optional parameters to send in the query
        :return:
        """
        try:
            if params:
                cursor.execute(query, params)
            else:
                cursor.execute(query)
        except psycopg2.Error as e:
            raise ValueError("Error executing SQL query: {}".format(e))
