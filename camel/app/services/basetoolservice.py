from collections import OrderedDict
from typing import Optional, List

import abc

from camel.app.connection.connection import Connection
from camel.app.parameter.parameter import Parameter
from camel.app.services.service import Service


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
    def get_dependencies(self) -> List[str]:
        """
        Gets the dependencies for the tool from the database
        :return: List with dependencies as string
        """
        pass

    @abc.abstractmethod
    def get_default_parameters(self) -> OrderedDict[str, Parameter]:
        """
        Returns the default parameters for this tool.
        :return: Default parameters
        """
        pass

    @abc.abstractmethod
    def get_names_mandatory_parameter(self) -> List[str]:
        """
        Returns all the default parameters.
        :return: Default parameter names
        """
        pass

    @abc.abstractmethod
    def get_parameter(self, parameter_name: str) -> Optional[Parameter]:
        """
        Returns the parameter with the given name.
        :param parameter_name: Parameter name
        :return: Parameter
        """
        pass
