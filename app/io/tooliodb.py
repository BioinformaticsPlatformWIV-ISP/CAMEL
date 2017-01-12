import os.path

from app.components.files.fileutils import FileUtils
from app.connection.connection import Connection
from app.io.toolio import ToolIO
from config import DB_CONFIG


class ToolIODb(ToolIO):
    """
    Class that represents an input / output database of a tool.
    """
    TYPE_NAME = 'db'

    def __init__(self, name, version='latest', logged=True, config=None):
        """
        Initializes a tool input / output database.
        :param name: Name of the database
        :param version: Version of the database
        :param logged: If True, the output can be logged
        :param config: Configuration file for the database connection
        """
        super(ToolIODb, self).__init__(logged)
        self._config = DB_CONFIG if config is None else config
        self._name = name
        self._version = version
        self._path = self.__get_location()

    @property
    def path(self):
        """
        Returns the value.
        :return: Value
        """
        return self._path

    def __str__(self):
        """
        String representation
        :return: String representation
        """
        return str(self.path)

    def __repr__(self):
        """
        Internal representation
        :return: Internal representation representation
        """
        return 'ToolIODb("{}")'.format(self.path)

    def is_valid(self):
        """
        Checks if the tool input / output file is valid.
        :return: True if valid
        """
        if not self.exists:
            return False
        return True

    def is_dir(self):
        """
        Checks whether the database is a directory or a file
        :return: True when the database is a directory
        """
        return True if os.path.isdir(self.path) else False

    @property
    def exists(self):
        """
        Checks whether this file exists.
        :return: True if the file exists, False otherwise
        """
        return os.path.isfile(self.path) if not self.is_dir() else True

    @property
    def hash(self):
        """
        Returns the hash value.
        :return: Hash value
        """
        return FileUtils.hash_directory(self.path) if self.is_dir() else FileUtils.hash_file(self.path)

    def __get_location(self):
        """
        Returns the location of the database
        :return: Database location
        """
        if self._version.lower() == 'latest':
            return self.__get_latest_loc()
        else:
            return self.__get_version_loc()

    def __get_latest_loc(self):
        """
        Gets the location of the 'latest' database
        :return: Database location
        """
        db_conn = Connection(self._config, 'db_loc')
        sql = """SELECT location FROM databases.db_loc
                    WHERE db_loc.name = %s
                      AND db_loc.latest IS TRUE
                      AND db_loc.available IS TRUE"""
        return self.__check_and_return_loc(db_conn.query(sql, (self._name, )))

    def __get_version_loc(self):
        """
        Gets the location of the database with the specified version
        :return: Database location
        """
        db_conn = Connection(self._config, 'db_loc')
        sql = """SELECT location FROM databases.db_loc
                    WHERE db_loc.name = %s
                      AND db_loc.version = %s
                      AND db_loc.available IS TRUE"""
        return self.__check_and_return_loc(db_conn.query(sql, (self._name, self._version)))

    def __check_and_return_loc(self, results):
        """
        Checks whether a single database is returned from the database. Otherwise an exception is thrown.
        :param results: Array with the results of the database query
        :return: Database location
        """
        if len(results) == 1:
            raise ValueError('No available database found with name {} and version {}'.format(self._name,
                                                                                              self._version))
        elif len(results) > 2:
            raise ValueError('Too many results found for database with name {} and version {}'.format(self._name,
                                                                                                      self._version))
        else:
            return results[1][0]
