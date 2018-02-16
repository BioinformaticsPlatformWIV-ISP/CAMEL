from collections import OrderedDict
from typing import List, Optional

from app.connection.connection import Connection
from app.parameter.parameter import Parameter
from app.services.basetoolservice import BaseToolService


class ToolService(BaseToolService):
    """
    This class will perform operations on the DB regarding tools (retrieve tool parameters, info, etc).
    """

    def __init__(self, tool_name: str, tool_version: str, connection: Connection) -> None:
        """
        Initializes tool service class assigning the module
        :param tool_name: Name of the tool
        :param tool_version: Version of the tool
        :param connection: Connection to the database
        :return: none
        """
        super(ToolService, self).__init__(connection)
        self._tool_id = self.__get_tool_id(tool_name, tool_version)

    @property
    def tool_id(self) -> int:
        """
        Returns the id of the tool in the database.
        :return: Id
        """
        return self._tool_id

    def __get_tool_id(self, name: str, version: str) -> int:
        """
        Returns the tool id for this tool.
        :return: Tool id
        """
        sql = """
        SELECT tool_id FROM tools.tool
        WHERE LOWER(tool_func_name) = %s
        AND version = %s"""
        try:
            return self.db_connection.query(sql, (name.lower(), version))[1][0]
        except IndexError:
            raise ValueError(f"Tool {name} {version} not found in the database.")

    def get_tool_command(self) -> str:
        """
        Returns the tool command specified in the database
        : return: Tool command
        """
        sql = """
        SELECT command FROM tools.tool
        WHERE tool_id = %s;
        """
        return self.db_connection.query(sql, (self._tool_id,))[1][0]

    def get_dependencies(self) -> Optional[str]:
        """
        Gets the dependencies for the tool from the database
        : return: String with whitespace separated list of dependencies
        """
        sql = """
        SELECT d.lmod
        FROM tools.tool_dependency td, tools.dependency d
        WHERE tool_id = %s
        AND td.dependency_id = d.dependency_id;
        """
        dependencies = self.db_connection.query(sql, (self._tool_id,))
        return None if len(dependencies) == 1 else ' '.join([dep[0] for dep in dependencies[1:]])

    def get_default_parameters(self) -> OrderedDict:
        """
        Returns the default parameters for this tool.
        :return: Default parameters
        """
        sql = """
        SELECT name, option, value, COALESCE(p_index, 0) as p_index
        FROM tools.tool_parameter
        WHERE tool_id = %s
        AND active = TRUE;
        """
        query_result = self.db_connection.query(sql, (self._tool_id,))
        parameters = OrderedDict()
        for name, option, value, _ in sorted(query_result[1:], key=lambda x: x[3]):
            parameters[name] = Parameter(name, option, value)
        return parameters

    def get_all_parameters(self) -> OrderedDict:
        """
        Returns all parameters for this tool.
        :return: All parameters
        """
        sql = """
        SELECT name, option, value, COALESCE(p_index, 0) as p_index
        FROM tools.tool_parameter
        WHERE tool_id = %s;
        """
        query_result = self.db_connection.query(sql, [self._tool_id])
        parameters = OrderedDict()
        for name, option, value, _ in sorted(query_result[1:], key=lambda x: x[3]):
            parameters[name] = Parameter(name, option, value)
        return parameters

    def get_names_mandatory_parameter(self) -> List[str]:
        """
        Returns all the default parameters.
        :return: Default parameter names
        """
        sql = """
        SELECT name FROM tools.tool_parameter
        WHERE tool_id = %s
        AND mandatory = True;"""
        query_result = self.db_connection.query(sql, (self._tool_id,))
        return [x[0] for x in query_result[1:]]

    def get_parameter(self, parameter_name) -> Optional[Parameter]:
        """
        Returns the parameter with the given name.
        :param parameter_name: Parameter name
        :return: Parameter
        """
        sql = """
        SELECT name, option, value, p_index
        FROM tools.tool_parameter
        WHERE tool_id = %s
        AND name = %s;
        """
        try:
            parameter_row = self.db_connection.query(sql, (self._tool_id, parameter_name))[1]
            return Parameter(parameter_row[0], parameter_row[1], parameter_row[2])
        except IndexError:
            return None
