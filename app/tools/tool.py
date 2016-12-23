import abc
import copy

import logging
import os

from app.command.command import Command
from app.io.tool_io import ToolIO
from app.loggers.log_manager import LogManager
from app.services.tool_service import ToolService


class Tool(object):
    """
    Contains the common functionality of tools.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, name, version, camel):
        """
        Initializes a tool.
        """
        self._name = name
        self._version = version
        self._tool_inputs = {}
        self._tool_outputs = {}
        self._tool_service = ToolService(name, version, camel.connection)
        self._tool_command = self._tool_service.get_tool_command()
        self._parameters = self._tool_service.get_default_parameters()
        self._command = Command()
        self._folder = None
        self._log_handler = None

    @property
    def name(self):
        """
        Returns the name of this tool.
        :return: Name
        """
        return '{} {}'.format(self._name, self._version)

    @property
    def stdout(self):
        """
        Returns the command line stdout.
        :return: Stdout
        """
        return self._command.stdout

    @property
    def stderr(self):
        """
        Returns the command line stderr.
        :return: Stderr
        """
        return self._command.stderr

    def add_input_files(self, input_files):
        """
        Sets the input files for a tool.
        :param input_files: New input files
        :return: None
        """
        self._tool_inputs.update(input_files)

    def update_parameters(self, **kwargs):
        """
        Updates the parameters for this tool.
        :param kwargs: Parameters in key value format
        :return: None
        """
        for parameter_name, new_value in kwargs.iteritems():
            parameter = self._tool_service.get_parameter(parameter_name)
            if not parameter:
                raise ValueError("{} has no parameter '{}'".format(self._name, parameter_name))
            if new_value is False:
                del(self._parameters[parameter_name])
            else:
                parameter.value = str(new_value)
                self._parameters[parameter_name] = parameter

    def run(self, folder='.'):
        """
        Runs this tool.
        :param folder: Folder to run the tool in.
        :return: None
        """
        self._folder = folder
        self._log_handler = LogManager.get_file_handler(os.path.join(self._folder, 'debug.log'), logging.DEBUG)
        logging.getLogger().addHandler(self._log_handler)
        logging.info('Running tool {}'.format(self.name))
        self._check_parameters()
        self._check_input()
        self._execute_tool()
        self._check_output()
        logging.getLogger().removeHandler(self._log_handler)

    def get_outputs(self, key):
        """
        Returns the outputs with the given key.
        :param key: output key
        :return: Output list
        """
        if key not in self._tool_outputs:
            raise ValueError("No output file with key '{}' found".format(key))
        return copy.copy(self._tool_outputs[key])

    def _build_dependencies(self):
        """
        Builds the dependencies.
        :return: Command to load dependencies
        """
        dependencies = self._tool_service.get_dependencies()
        if dependencies is not None:
            return 'module load ' + dependencies + '; '
        else:
            return ''

    def _build_options(self, excluded_parameters=None, delimiter=' '):
        """
        Builds the options string.
        :parameter delimiter: Delimiter between option and value
        :return: Options string
        """
        options = []
        for name, parameter in self._parameters.iteritems():
            if (excluded_parameters is not None) and (name in excluded_parameters):
                continue
            if parameter.value is not None:
                options.append(parameter.option + delimiter + str(parameter.value))
            else:
                options.append(parameter.option)
        return options

    def _execute_command(self):
        """
        Executes a the command.
        :return: None
        """
        if self._command.command is None:
            raise ValueError("Command is 'None'.")
        self._command.command = self._build_dependencies() + self._command.command
        self._command.run_command(self._folder)

    @abc.abstractmethod
    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        pass

    def _check_parameters(self):
        """
        Checks if the tool parameters are valid.
        :return: None
        """
        mandatory_parameters = self._tool_service.get_names_mandatory_parameter()
        for mandatory_parameter in mandatory_parameters:
            if mandatory_parameter not in self._parameters:
                raise ValueError("Mandatory parameter {} not set".format(mandatory_parameter))

    def _check_input(self):
        """
        Checks if the tool input is valid.
        :return: None
        """
        for output_key, output_list in self._tool_inputs.iteritems():
            for tool_output in output_list:
                if tool_output is None:
                    raise ValueError("Tool input with key {} is None".format(output_key))
                if not tool_output.is_valid():
                    raise ValueError("Invalid tool input with key {}: {}".format(output_key, tool_output))

    def _check_output(self):
        """
        Checks if the output is valid.
        :return: None
        """
        for output_key, output_list in self._tool_outputs.iteritems():
            for tool_output in output_list:
                if tool_output is None:
                    raise ValueError("Tool output with key {} is None".format(output_key))
                if not isinstance(tool_output, ToolIO):
                    raise ValueError("'{} {}' is not a tool output object".format(tool_output, type(tool_output)))
                if not tool_output.is_valid():
                    raise ValueError("Invalid tool output with key {}: {}".format(output_key, tool_output))

    def display(self):
        """
        Displays the current state of the tool.
        :return: None
        """
        print('{} {}'.format(self._name, self._version))
        print('-'*20)
        print('- Input -')
        for key, value in self._tool_inputs.iteritems():
            print('{}: {}'.format(key, ', '.join([str(o) for o in value])))
        print('')
        print('- Parameters -')
        for name, parameter in self._parameters.iteritems():
            print('{}: {}'.format(name, parameter))
        print('')
        print('- Output -')
        for key, value in self._tool_outputs.iteritems():
            print('{}: {}'.format(key, ', '.join([str(o) for o in value])))
        print('')
        print('-'*20)
        print('')
