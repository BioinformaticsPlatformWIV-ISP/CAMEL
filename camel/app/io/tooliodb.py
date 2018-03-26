import os.path

from camel.app.components.files.fileutils import FileUtils
from camel.app.connection.connection import Connection
from camel.app.io.toolio import ToolIO
from camel.config import DB_CONFIG


class ToolIODb(ToolIO):
    """
    Class that represents an input / output database of a tool.
    """

    def __init__(self, name: str, version: str='latest', logged: bool=True, config: str=None, prefix: str=None):
        """
        Initializes a tool input / output database.
        :param name: Name of the database
        :param version: Version of the database
        :param logged: If True, the output can be logged
        :param config: Configuration file for the database connection
        :param prefix: Prefix for the database to append to the path
        """
        super(ToolIODb, self).__init__(logged)
        self._config = DB_CONFIG if config is None else config
        self._name = name
        self._version = version
        self._prefix = prefix
        self._path = self.__get_location()

    @property
    def path(self) -> str:
        """
        Returns the value.
        :return: Value
        """
        return self._path

    def __str__(self) -> str:
        """
        String representation
        :return: String representation
        """
        return str(self.path)

    def __repr__(self) -> str:
        """
        Internal representation
        :return: Internal representation representation
        """
        return f'ToolIODb("{self.path}")'

    def is_valid(self) -> bool:
        """
        Checks if the tool input / output file is valid.
        :return: True if valid
        """
        return self.exists

    def is_dir(self) -> bool:
        """
        Checks whether the database is a directory or a file
        :return: True when the database is a directory
        """
        return True if os.path.isdir(self.path) else False

    @property
    def exists(self) -> bool:
        """
        Checks whether this file exists.
        :return: True if the file exists, False otherwise
        """
        if self._prefix is None:
            return True if self.is_dir() else os.path.isfile(self.path)
        else:
            return os.path.isdir(os.path.dirname(self._path))

    @property
    def hash(self) -> str:
        """
        Returns the hash value.
        :return: Hash value
        """
        return FileUtils.hash_directory(self.path) if self.is_dir() else FileUtils.hash_file(self.path)

    @property
    def type_name(self) -> str:
        """
        Returns the type of the IO object.
        :return: Type value
        """
        return 'db'

    def __get_location(self) -> str:
        """
        Returns the location of the database
        :return: Database location
        """
        return self.__get_db_location(self._version.lower()) if self._prefix is None else os.path.join(self.__get_db_location(self._version.lower()), self._prefix)

    def __get_db_location(self, version: str, prefix: str=None) -> str:
        """
        Gets the location of the database in the file system
        :param version: Version of the database
        :param prefix: Prefix for the database
        :return: Database location in the file system
        """
        db_conn = Connection(self._config, 'db_loc')
        sql = """SELECT location FROM databases.db_loc 
                    WHERE db_loc.name = %s 
                    AND db_loc.available IS TRUE"""
        sql += " AND db_loc.latest IS TRUE" if (version == 'latest') else " AND db_loc.version = %s"
        if prefix is None:
            return self.__check_and_return_loc(db_conn.query(sql, (self._name,)))
        else:
            return self.__check_and_return_loc(db_conn.query(sql, (self._name, self._version)))

    def __check_and_return_loc(self, results) -> str:
        """
        Checks whether a single database is returned from the database. Otherwise an exception is thrown.
        :param results: Array with the results of the database query
        :return: Database location
        """
        if len(results) != 2:
            raise ValueError(f'There can be only one database with name {self._name} and version {self._version}, but {len(results)-1} were found!')
        else:
            return results[1][0]
