from collections import OrderedDict
from typing import Optional

import abc

from app.connection.connection import Connection
from app.parameter.parameter import Parameter
from app.services.service import Service


class BaseToolService(Service, metaclass=abc.ABCMeta):
    """
    This is the base class for the tool service.
    """

    def __init__(self, connection: Optional[Connection]) -> None:
        """
        Initializes the base tool service.
        :param connection: Connection to the database
        :return: None
        """
        super().__init__(connection)

    @property
    @abc.abstractmethod
    def tool_id(self) -> int:
        """
        Returns the id of the tool in the database.
        :return: Id
        """
        pass

    @abc.abstractmethod
    def get_tool_command(self) -> str:
        """
        Returns the tool command specified in the database
        :return: Tool command
        """
        pass

    @abc.abstractmethod
    def get_dependencies(self) -> Optional[str]:
        """
        Gets the dependencies for the tool from the database
        :return: String with whitespace separated list of dependencies
        """
        pass

    @abc.abstractmethod
    def get_default_parameters(self) -> OrderedDict:
        """
        Returns the default parameters for this tool.
        :return: Default parameters
        """
        pass

    # def get_all_parameters(self) -> OrderedDict:
    #     """
    #     Returns all parameters for this tool.
    #     :return: All parameters
    #     """
    #     sql = """
    #     SELECT name, option, value, COALESCE(p_index, 0) as p_index
    #     FROM tools.tool_parameter
    #     WHERE tool_id = %s;
    #     """
    #     query_result = self.db_connection.query(sql, [self._tool_id])
    #     parameters = OrderedDict()
    #     for name, option, value, _ in sorted(query_result[1:], key=lambda x: x[3]):
    #         parameters[name] = Parameter(name, option, value)
    #     return parameters

    # def get_names_mandatory_parameter(self) -> List[str]:
    #     """
    #     Returns all the default parameters.
    #     :return: Default parameter names
    #     """
    #     sql = """
    #     SELECT name FROM tools.tool_parameter
    #     WHERE tool_id = %s
    #     AND mandatory = True;"""
    #     query_result = self.db_connection.query(sql, (self._tool_id,))
    #     return [x[0] for x in query_result[1:]]

    @abc.abstractmethod
    def get_parameter(self, parameter_name) -> Optional[Parameter]:
        """
        Returns the parameter with the given name.
        :param parameter_name: Parameter name
        :return: Parameter
        """
        pass
